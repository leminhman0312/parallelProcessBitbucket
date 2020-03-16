#include <cstdlib>
#include <fstream>
#include <mpi.h>
#include <cmath>
#include <iostream>

double *u, *x, *F, *us;
double a = 0.1;
double t, dt, tf;

#define NGHOST para.nghost
#define NX 21
#define XMIN 0.
#define XMAX 1.
#define DX ((XMAX-XMIN)/(NX-1))

#define U(i) u[i]
#define ILO  para.ilo
#define IHI  para.ihi

#define ILOPAD 0
#define IHIPAD IHI+NGHOST

/* MPI STUFF */
struct mpistr {
    int ilo, ihi, np, myid;
    int idim, idip;
    int xs, xe, nghost;
    MPI_Comm comm;
};

mpistr para;

double *sendleft, *sendright, *recleft, *recright;
MPI_Request sendrl, sendrr;
MPI_Request recrl, recrr;

#define NP      para.np
#define IAMNODE para.myid
#define NRIGHT  para.idip
#define NLEFT   para.idim


/* methods */
void boundary();
void setup_MPI();
void update(int *flag);
void bc_send();
void bc_receive();
void sync_nodes();
void saveout(int ii);

int main(int argc, char* argv[]) {

    para.nghost = strtol(argv[1],&argv[1],10);

    setup_MPI();
    
    MPI_Barrier(para.comm);

    if(IAMNODE==0) fprintf(stdout, "Running with %3d procs and %3d ghost\n", NP, NGHOST);

    double xc = 0.05; /* width of the gaussian */

    /* initialize */
    u = new double[IHIPAD]();
    F = new double[IHIPAD+1]();
    x = new double[NX+2*NGHOST]();

    us = new double[NX]();

    sendleft  = new double[NGHOST]();
    sendright = new double[NGHOST]();
    recleft   = new double[NGHOST]();
    recright  = new double[NGHOST]();



    for(int i=0; i<NX+2*NGHOST; i++) {
        x[i] = XMIN + ((double) (i-NGHOST)*DX);
    }

    /* initial condition */
    for(int i=ILOPAD; i<IHIPAD; i++) {
        int ii = i + para.xs;
        double xi = x[ii];

        U(i) = exp(-(xi-0.5)*(xi-0.5)/(2.*xc*xc));

    }

    /* runtime parameters */
    dt = 1.e-2;
    t  = 0.;
    tf = 2.;

    int ii=0;
    int savecounter = 0;
    int flag = 0;

    while(t < tf) {
        update(&flag);
        printf("This is flag : %d",flag);
        // flag = 0;
        t += dt;
        // fprintf(stdout,"%.3f\n", t);

        if((ii%5)==0) {
            saveout(savecounter++);
        }

        ii++;

    }

    if(IAMNODE==0) fprintf(stdout,"Ran %3d steps\n", ii);

    // update();

    MPI_Finalize();
    return 0;
}

void boundary() {
    if(IAMNODE==0) {
        U(1) = U(2);
        U(0) = U(1);
    }

    if(IAMNODE==(NP-1)) {
        U(IHI)   = U(IHI-1);
        U(IHI+1) = U(IHI);
    }
}


void update(int *flag) {

    int npoint;
    if(*flag == 0) {
        bc_send();
        bc_receive();

        *flag = NGHOST;
    }
    
    npoint = *flag-1;

    for(int i=ILO-npoint-1; i<IHI+npoint+1; i++) {
        if( a >= 0.) {
            F[i] = a*U(i-1);
        } else {
            F[i] = a*U(i);
        }
    }


    for(int i=ILO-npoint; i<IHI+npoint; i++) {
        U(i) -= dt*(F[i+1] - F[i])/DX;
    }

    boundary();

    *flag -= 1;

}

void setup_MPI() {
    // Start MPI
    MPI_Init(NULL, NULL);
    para.comm = MPI_COMM_WORLD;
    // Get the number of processors
    MPI_Comm_size(para.comm, &para.np);
    // Get the current processor
    MPI_Comm_rank(para.comm, &para.myid);

    // Figure out where the node configuration is
    int nblocks, nrem;
    nblocks = floor(NX/NP);
    nrem    = NX - nblocks*NP;
    
    // If it's the last node, add the 
    // add the remaining blocks 
    ILO = NGHOST;
    IHI = nblocks + NGHOST;
    para.xs = nblocks*IAMNODE;
    para.xe = nblocks*(IAMNODE+1);

    if(IAMNODE==(NP-1)) {
        IHI += nrem;
        para.xe += nrem;
    }



    // Figure out the processor before
    NLEFT = IAMNODE - 1;
    // Figure out the processor after
    NRIGHT = IAMNODE + 1;

    fprintf(stdout,"I am %d and ILO=%3d and IHI=%3d and xs=%3d and xe=%3d\n", IAMNODE, ILO, IHI, para.xs, para.xe);
}


void bc_send() {

    if(IAMNODE>0) {
        MPI_Irecv(recleft, NGHOST, MPI_DOUBLE, NLEFT, 1, para.comm, &recrl);
    }

    if(IAMNODE<(NP-1)) {
        MPI_Irecv(recright, NGHOST, MPI_DOUBLE, NRIGHT, 2, para.comm, &recrr);
    }

    for(int j=0; j<NGHOST; j++) {
        sendleft[j]  = U(ILO+j);
        sendright[j] = U(IHI-j-1);
    }

    if(IAMNODE > 0) {
        MPI_Isend(sendleft, NGHOST, MPI_DOUBLE, NLEFT, 2, para.comm, &sendrl);
    }

    if(IAMNODE<(NP-1)) {
        MPI_Isend(sendright, NGHOST, MPI_DOUBLE, NRIGHT, 1, para.comm, &sendrr);
    }
}


void bc_receive() {
    if(IAMNODE>0) {
        MPI_Wait(&sendrl, MPI_STATUS_IGNORE);
        MPI_Wait(&recrl, MPI_STATUS_IGNORE);
    }

    if(IAMNODE<(NP-1)) {
        MPI_Wait(&sendrr, MPI_STATUS_IGNORE);
        MPI_Wait(&recrr, MPI_STATUS_IGNORE);
    }

    for(int j=0; j<NGHOST; j++) {
        U(ILO-j-1) = recleft[j];
        U(IHI+j)   = recright[j];
    }
}

void sync_nodes() {
    double *unode;

    if(IAMNODE==0) {
        for(int p=1; p<NP; p++) {
            int ilo, ihi, nx;
            

            MPI_Recv(&ilo, 1, MPI_INT, p, 0, para.comm, MPI_STATUS_IGNORE);
            MPI_Recv(&ihi, 1, MPI_INT, p, 1, para.comm, MPI_STATUS_IGNORE);

            nx = ihi - ilo;

            unode = new double[nx]();

            MPI_Recv(unode, nx, MPI_DOUBLE, p, 2, para.comm, MPI_STATUS_IGNORE);


            for(int i=ilo; i<ihi; i++) {
                int ii = i - ilo;
                us[i] = unode[ii];
            }
            
        }

        for(int i=ILO; i<IHI; i++) {
            int ii = i-ILO;
            us[ii] = U(i);
        }
    } else {
        MPI_Send(&para.xs, 1, MPI_INT, 0, 0, para.comm);
        MPI_Send(&para.xe, 1, MPI_INT, 0, 1, para.comm);

        unode = new double[IHI-ILO]();

        for(int i=ILO; i<IHI; i++) {
            int ii = i-ILO;
            unode[ii] = U(i);
        }

        MPI_Send(unode, IHI-ILO, MPI_DOUBLE, 0, 2, para.comm);
    }

}


void saveout(int ii) {

    sync_nodes();

    
    if(IAMNODE==0) {
        char filename[50];

        sprintf(filename,"save/spook%03d.csv", ii);
        FILE * outfile; 

        outfile = fopen(filename, "w");

        for(int i=0; i<NX; i++) {
            fprintf(outfile, "%.3f,%.3f\n",x[i+NGHOST], us[i]);
        }
        fclose(outfile);
    }





}
