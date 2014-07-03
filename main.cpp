#include <iostream>
#include <ctime>
#include <queue>
#include "config.h"
#include "PartArray.h"
#include "mpi.h"

using namespace std;

void moveSystemPosRandomly(PartArray* sys, double d){
    Vect dir; //направление, в которое двигать частицу.
    for(int i=0;i<sys->count();i++){
        double longitude = ((double)config::Instance()->rand()/(double)config::Instance()->rand_max) * 2. * M_PI;
        dir.x = d * cos(longitude);
        dir.y = d * sin(longitude);
        dir.z = 0;
        sys->parts[i]->pos += dir;
    }
}

void moveSystemMRandomly(PartArray* sys, double fi){
    double side = 1.;
    double oldFi;

    for(int i=0;i<sys->count();i++){
        double oldLen = sys->parts[i]->m.length();
        if ((double)config::Instance()->rand()/(double)config::Instance()->rand_max>0.5)
            side = -1.;
        else
            side = 1.;
        oldFi = sys->parts[i]->m.angle();
        double longitude = oldFi+(fi*side);
        sys->parts[i]->m.x = oldLen * cos(longitude);
        sys->parts[i]->m.y = oldLen * sin(longitude);
        sys->parts[i]->m.z = 0;
    }
}

int main(int argc, char** argv)
{
    MPI_Init (&argc, &argv);	/* starts MPI */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size==1){
        cout<<"invalid process count"<<endl;
        MPI_Finalize();
        return 0;
    }

    config::Instance()->set3D();
    config::Instance()->srand(time(NULL)+rank*size);config::Instance()->rand();

    //директивы
    int x=3, y=3, z=3, //количество частиц в линейке
            experimentCount=10;
    double intervalCount = 100.;
    double space = config::Instance()->partR*4.;//расстояние между центрами частиц
    double dMax = space/2.-config::Instance()->partR;
    int exitCode = -1; //код выхода из системы
    double transBuff[4] = {0,0,0,0}; //буффер передачи. 0 - d, 1 - m, 2 - n(количество поворотов), 3 - с (количество попыток);

    PartArray *sys, *example;
    example = new PartArray((double)x*space,(double)y*space,(double)z*space);//размер образца 4 радиуса, или 2 диаметра

    //бросаем частицы в шахматном порядке на линию и запоминаем состояние
    example->dropChain(space);
    StateMachineFree oldState(example->state);

    if (rank!=0){
        //все задачи берут свои рассчеты и выполняют
        int iteration = 0; //номер итерации
        double d=0;
        for (int i=0; i<intervalCount; i++){
            double m=0.;
            for (int k=0;k<intervalCount;k++){
                if (iteration==rank-1){


                    //сам эксперимент
                    int anomCount=0;
                    for (int j=0;j<experimentCount;j++){
                        sys = example->copy();

                        moveSystemMRandomly(sys,m);
                        moveSystemPosRandomly(sys,d);

                        sys->setToGroundState();

                        if (!(oldState==sys->state)){
                            anomCount++;
                        }

                        delete sys;
                    }

                    //отправляем данные в 1 поток
                    transBuff[0] = d;
                    transBuff[1] = m;
                    transBuff[2] = (double)anomCount;
                    transBuff[3] = (double)experimentCount;
                    MPI_Send(&transBuff,4,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
                    //cout<<rank<<": send code "<<anomCount<<endl;

                }
                if (iteration==size-2)
                    iteration=0;
                else
                    iteration++;

                m+=M_PI_2/intervalCount;
            }
            d+=dMax/intervalCount;
        }


        //отправляем код завершения процесса
        transBuff[0] = exitCode;
        transBuff[1] = transBuff[2] = transBuff[3] = 0;
        MPI_Send(&transBuff,4,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        //cout<<rank<<": send code "<<exitCode<<endl;

    } else {
        ofstream f("checkMagnetism2DRes_3x7.dat");
        f<<"d\tfi\tanom\tcount"<<endl;

        MPI_Status stat;
        int finished=0, transactionCount=0;
        do{
            MPI_Recv(&transBuff,4,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&stat);
            transactionCount++;
            cout << transactionCount/(intervalCount*intervalCount+size-1)*100<<"% complete"<<endl;
            //cout<<rank<<": recieve code "<<transBuff[2]<<" from "<<stat.MPI_SOURCE<<endl;
            if (transBuff[0]==exitCode)
                finished++;
            else {

                //сохраняем данные
                f<<
                    transBuff[0]<<"\t"<<
                      transBuff[1]<<"\t"<<
                      transBuff[2]<<"\t"<<
                      transBuff[3]<<"\t"<<endl;

            }
        } while (finished<(size-1));
        cout<<rank<<": "<<"program may finish"<<endl;
        f.close();
    }

    MPI_Finalize();
    return 0;
}

