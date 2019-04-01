#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.14159265358979323846

float contForces2(float cc_part_i[3][3], float U_part_i[4][4], float V_part_i[3][3],  float MatPart[3][3],
                  float UPart[3][3], float VPart[3][3], float NormalF[4], float TangF[6], float ks,
                  float kn, float ro, float dt, int partIDs[6]){
                      double mmnn,ll,NormalF_p1[11][1],TangF_p1[11][1],v_Forces[3][sizeof(mmnn)];
                      mmnn = sizeof(MatPart);
                      ll = sizeof(NormalF);
                      int i,j,nnz[9],ct;
                      ct = 0;
                      for (i=0; i<3; ++i){
                          for (j=0;j<3;++j){
                                if (UPart[i][j]!=0){
                                    UPart[i][j] = nnz[ct];
                                    ct+=1;
                                    return nnz[9];
                                }
                          }
                      }

                      int k,p,nnzz[16],cont;
                      cont = 0;
                      for (k=0; k<4; ++k){
                          for (p=0;p<4;++p){
                                if (U_part_i[k][p]!=0){
                                    U_part_i[k][p] = nnzz[cont];
                                    cont+=1;
                                    return nnzz[16];
                                }
                          }
                      }
                      if (nnz[9]==0 && nnzz[16]==0){
                        NormalF_p1[11][1] = NormalF[4];
                        TangF_p1[11][1] = TangF[6];
                        return 0;
                      }

                      int t,o,q,w;
                      float Fs_ijmax,Fn_ij,Fs_ij,dfs,dfn,vs,vn,e[sizeof(mmnn)],d2,trans[3][3],rest[3][3],Fnt[6],Fst[6],A[6];
                      double  cc_part_j[sizeof(MatPart)][sizeof(MatPart)],U_part_j[sizeof(MatPart)][sizeof(MatPart)],V_part_j[sizeof(MatPart)][sizeof(MatPart)];
                      rest[3][3] = (cc_part_j[sizeof(MatPart)][sizeof(MatPart)]-cc_part_i[3][3]);
                      for (q = 0; q<mmnn;++q){
                        for (w = 0; w<mmnn; ++w){
                            trans[w][q]= rest [q][w];
                        }
                      }
                      for (t=0; t < mmnn; ++t){
                          for (o=0; o<mmnn; ++o){
                                cc_part_j[sizeof(MatPart)][sizeof(MatPart)] = MatPart[t][o];
                                U_part_j[sizeof(MatPart)][sizeof(MatPart)] = UPart[t][o];
                                V_part_j[sizeof(MatPart)][sizeof(MatPart)] = VPart[t][o];

                          if( (cc_part_j[1][1]==cc_part_i[1][1]) && (cc_part_j[2][2]==cc_part_i[2][2])){

                          }
                            else{
                                d2 = (cc_part_j[sizeof(MatPart)][sizeof(MatPart)]-cc_part_i[3][3])*(trans[3][3]);

                                Fnt[6] = NormalF[partIDs[t]];
                                Fst[6] = TangF[partIDs[t]];
                                if (sqrt(d2)<= 2 * ro) {

                                    vn = e[1] * (V_part_i[1][1] - V_part_j[1][1]) + e[2] * (V_part_i[2][2] - V_part_j[2][2]);
                                    vs = e[2] * (V_part_i[1][1] - V_part_j[1][1]) - e[2] * (V_part_i[2][2] - V_part_j[2][2]) - ro * (V_part_i[3][3] + V_part_j[3][3])/1000.0;

                                    dfn = kn * vn * dt;
                                    dfs = ks * vs * dt;

                                    Fn_ij = Fnt[6] + dfn;
                                    Fs_ij = Fst[6] + dfs;

                                    if (Fn_ij < 0){
                                        Fn_ij = 0;
                                        Fs_ij = 0;
                                    }
                                    Fs_ijmax = Fn_ij * tan((30.0*M_PI)/180.0);

//                                    if abs(Fs_ij > Fs_ijmax){
//                                        Fs_ij = abs(Fs_ijmax)*Fs_ij; /*Falta el sign en C*/
                                    }
//                                    NormalF_p1[partIDs[i]] = Fn_ij;
//                                    TangF_p1[partIDs[i]]  =  Fs_ij;




                                }
                          }
                          }
                      }












////                      return v_Forces, NormalF_p1, TangF_p1;

/*Every integer has a weight of 4*/
//
////    int i, j, count = 0;
////    for (i = 0; i <  r; i++)
////      for (j = 0; j < c; j++)
////         *(arr + i*c + j) = ++count;
//    int i, j;
//    for (i = 0; i <  r; i++)
////        for (j = 0; j < c; j++)
//            printf("%d ", arr[i]);
//
//   /* Code for further processing and free the
//      dynamically allocated memory */
//
//   return 0;
//}
//    for(int i = 0; i < 2; i++) {
//        for(int j = 0; j < 2; j++) {
//            printf("%d ", arr[i][j]);
//    }
//    printf("\n");



//
//
//    int row =sizeof(MatPart)/sizeof(MatPart[0]);
//    int column = sizeof(MatPart[0])/sizeof(MatPart[0][0]);
//    printf("%d\n", sizeof(MatPart));
//    printf("%d\n", sizeof(MatPart[0][0]));
//
//
//    printf("Number of rows: %d\n", row);
//    printf("Number of columns: %d\n", column);




//int column = sizeof(MatPart[0])/sizeof(MatPart[0][0]);

//
//                    float mmnn sizeof(MatPart);
//                    float ll [NormalF][1];
//                    float NormalF_p1[11][1];
//                    float TangF_p1 [11][1];
//                    float v_Forces [3][MatPart];
//
//                    int Num;
//                    int i;
//                    for (i=0; i<MatPart; ++i){
//                        if (UPart != 0 ){
//                            Num = Num + 1;
//                        }
//                    }
//                    int j;
//                    public static double cc_part_k[][] = new double[MatPart][MatPart];
//                    double U_part_k;
//                    double V_part_k;
//                    for (j = 0 ; j < MatPart; ++j){
//                        cc_part_k = MatPart[j][j];
//                        double U_part_k[] = MatPart[j][j];
//                        double V_part_k[] = MatPart[j][j];
//
//#define
//                    }

//                    for i=1:mm
//
//                    cc_part_j = MatPart(i,:);
//                    U_part_j  = UPart(i,:);
//                    V_part_j  = VPart(i,:);
//
//                    int j;
//                    double cc_part_k[MatPart][MatPart];
//                    double U_part_k[MatPart][MatPart];
//                    double V_part_k[MatPart][MatPart];
//                    for (j; j<MatPart; ++j){
//                        cc_part_k[MatPart][MatPart]= MatPart[j][j];
//
//                        if((cc_part_k[1]==cc_part_i[1]) && (cc_part_k[2]==cc_part_i[2]))
//                            v_Forces[][] = 0;
//                        else
//                            d2 = (cc_part_k-cc_part_i)*(cc_part_k-cc_part_i);
//                        Fnt = NormalF[partIDs[i]];
//                        Fst = Tang[partIDs[i]];
//
//                        if sqrt(d2)<=2*ro
//                            e = (cc_part_k-cc_part_i)
//                    }
//
//                  }
