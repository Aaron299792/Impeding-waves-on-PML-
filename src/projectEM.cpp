/*******************************************************************
    Universidad de Costa Rica
    Projecto Electromagnetismo II -- IIS-2025
    Title: Impeding Electromagnetic Waves on a PML.
    Author: Aaron Sanabria Martínez. C17208
    
    Credits: This work was bases on the Python solution presented
    by Jennifer E. Houle, Dennis M. Sullivan (2019) in their book 
    Electromagnetic Simulation Using the FDTD Method with Python
*******************************************************************/

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fstream>
#include <string>

#ifdef _OPENMP
    #include <omp.h>
#endif

typedef std::vector<double> rvec;
typedef std::vector< std::complex<double> > cvec;

using namespace std;

constexpr double PI = 3.14159265358979323846;

void calculate_pml_parameters(int npml, int ie, int je, int ke, rvec& gi1, rvec& gi2, rvec& gi3, rvec& fi1, rvec& fi2, rvec& fi3, 
                              rvec& gj1, rvec& gj2, rvec& gj3, rvec& fj1, rvec& fj2, rvec& fj3, 
                              rvec& gk1, rvec& gk2, rvec& gk3, rvec& fk1, rvec& fk2, rvec& fk3){
    
    for (int n = 0; n < npml; ++n){
        double xxn = static_cast<double>(npml - n) / static_cast<double>(npml);
        double xn = 0.33 * pow(xxn, 3);
        fi1[n] = xn;
        fi1[ie - n - 1] = xn;
        gi2[n] = 1 / (1 + xn);
        gi2[ie - 1 - n] = 1 / (1 + xn);
        gi3[n] = (1 - xn) / (1 + xn);
        gi3[ie - 1 - n] = (1 - xn) / (1 + xn);

        fj1[n] = xn;
        fj1[je - n - 1] = xn;
        gj2[n] = 1 / (1 + xn);
        gj2[je - 1 - n] = 1 / (1 + xn);
        gj3[n] = (1 - xn) / (1 + xn);
        gj3[je - 1 - n] = (1 - xn) / (1 + xn);

        fk1[n] = xn;
        fk1[ke - n - 1] = xn;
        gk2[n] = 1 / (1 + xn);
        gk2[ke - 1 - n] = 1 / (1 + xn);
        gk3[n] = (1 - xn) / (1 + xn);
        gk3[ke - 1 - n] = (1 - xn) / (1 + xn);
        
        xxn = static_cast<double>(npml - n - 0.5) / static_cast<double>(npml);
        xn = 0.33 * pow(xxn, 3);

        gi1[n] = xn;
        gi1[ie - 1 - n] = xn;
        fi2[n] = 1 / (1 + xn);
        fi2[ie - n - 1] = 1 / (1 + xn);
        fi3[n] = (1 - xn) / (1 + xn);
        fi3[ie - n - 1] = (1 - xn) / (1 + xn);

        gj1[n] = xn;
        gj1[je - 1 - n] = xn;
        fj2[n] = 1 / (1 + xn);
        fj2[je - n - 1] = 1 / (1 + xn);
        fj3[n] = (1 - xn) / (1 + xn);
        fj3[je - n - 1] = (1 - xn) / (1 + xn);

        gk1[n] = xn;
        gk1[ke - 1 - n] = xn;
        fk2[n] = 1 / (1 + xn);
        fk2[ke - n - 1] = 1 / (1 + xn);
        fk3[n] = (1 - xn) / (1 + xn);
        fk3[ke - n - 1] = (1 - xn) / (1 + xn);
    }
}

void calculate_dx_fied(int ie, int je, int ke, rvec& dx, rvec& idx, rvec &hy, rvec &hz, rvec& gj3, rvec& gk3, rvec& gj2, rvec& gk2, rvec& gi1){
    #pragma omp for schedule(static) collapse(3)
    for (int i = 1; i < ie; ++i){
        for (int j = 1; j < je; ++j){
            for (int k = 1; k < ke; ++k){
                double curl_h =  hz[ (i*je + j)*ke + k ] -  hz[ ( i * je + (j - 1) )*ke + k ] - hy[ (i * je + j)*ke + k ] + hy[ (i*je + j)*ke + (k - 1) ];
                idx[ (i*je + j)*ke + k ] += curl_h;
                dx[ (i*je + j)*ke + k ] = gk3[k] * gj3[j] * dx[ (i*je + j)*ke + k ] + gj2[j] * gk2[k] * (0.5 * curl_h + gi1[i] * idx[ (i*je + j)*ke + k ]);
            }
        }
    }
}

void calculate_dy_fied(int ie, int je, int ke, rvec& dy, rvec& idy, rvec &hx, rvec &hz, rvec& gi3, rvec& gk3, rvec& gi2, rvec& gk2, rvec& gj1){
    #pragma omp for schedule(static) collapse(3)
    for (int i = 1 ; i < ie; ++i){
        for (int j = 1; j < je; ++j){
            for (int k = 1; k < ke; ++k){
                double curl_h =  hx[ (i*je + j)*ke + k ] - hx[ ( i*je + j )*ke + (k - 1) ] - hz[ (i*je + j)*ke + k ] + hz[ ((i - 1)*je + j)*ke + k ];
                idy[ (i*je + j)*ke + k ] += curl_h;
                dy[ (i*je + j)*ke + k ] = gi3[i] * gk3[k] * dy[ (i*je + j)*ke + k ] + gi2[i] * gk2[k] * (0.5 * curl_h + gj1[j] * idy[ (i*je + j)*ke + k ]); 
            }   
        }   
    }   
}

void calculate_dz_fied(int ie, int je, int ke, rvec& dz, rvec& idz, rvec &hx, rvec &hy, rvec& gi3, rvec& gj3, rvec& gi2, rvec& gj2, rvec& gk1){
    #pragma omp for schedule(static) collapse(3)
    for (int i = 1 ; i < ie; ++i){
        for (int j = 1; j < je; ++j){
            for (int k = 1; k < ke; ++k){
                double curl_h =  hy[ (i*je + j)*ke + k ] -  hy[ ( (i-1)*je + j )*ke + k ] - hx[ (i*je + j)*ke + k ] + hx[ (i*je + (j-1))*ke + k ];
                idz[ (i*je + j)*ke + k ] += curl_h;
                dz[ (i*je + j)*ke + k ] = gi3[i] * gj3[j] * dz[ (i*je + j)*ke + k ] + gi2[i] * gj2[j] * (0.5 * curl_h + gk1[k] * idz[ (i*je + j)*ke + k ]); 
            }   
        }   
    }   
}

void calculate_inc_dy_field(int ie, int je, int ke, int ia, int ib, int ja, int jb, int ka, int kb, rvec& dy, rvec& hx_inc){
    #pragma omp for schedule(static) collapse(2)
    for (int i = ia; i <= ib; ++i ){
        for (int j = ja; j <= jb; ++j){
            dy[ (i*je + j)*ke + ka ] -= 0.5*hx_inc[j];
            dy[ (i*je + j)*ke + kb + 1 ] += 0.5*hx_inc[j];
        }
    }
}

void calculate_inc_dz_field(int ie, int je, int ke, int ia, int ib, int ja, int jb, int ka, int kb, rvec& dz, rvec& hx_inc){
    #pragma omp for schedule(static) collapse(2)
    for (int i = ia; i <= ib; ++i ){
        for (int k = ka; k <= kb; ++k){
            dz[ (i*je + ja)*ke + k ] += 0.5*hx_inc[ja - 1];
            dz[ (i*je + jb)*ke + k ] -= 0.5*hx_inc[jb];
        }
    }
}

void calculate_e_fields(int ie, int je, int ke, rvec& dx, rvec& dy, rvec& dz, rvec& gax, rvec& gay, rvec& gaz, rvec& gbx, rvec& gby, 
        rvec& gbz, rvec& ex, rvec& ey, rvec& ez, rvec& ix, rvec& iy, rvec& iz){
 
    #pragma omp for schedule(static) collapse(3)
    for (int i = 0; i < ie; ++i){
        for (int j = 0; j < je; ++j){
            for (int k = 0; k < ke; ++k){
                ex[ (i*je + j)*ke + k ] = gax[ (i*je + j)*ke + k ]* (dx[ (i*je + j)*ke + k ] - ix[ (i*je + j)*ke + k ]);
                ix[ (i*je + j)*ke + k ] += gbx[ (i*je + j)*ke + k ]*ex[ (i*je + j)*ke + k ];
                ey[ (i*je + j)*ke + k ] = gay[ (i*je + j)*ke + k ]* (dy[ (i*je + j)*ke + k ] - iy[ (i*je + j)*ke + k ]);
                iy[ (i*je + j)*ke + k ] += gby[ (i*je + j)*ke + k ]*ey[ (i*je + j)*ke + k ];
                ez[(i*je + j)*ke + k ] = gaz[ (i*je + j)*ke + k ]* (dz[ (i*je + j)*ke + k ] - iz[ (i*je + j)*ke + k ]);
                iz[ (i*je + j)*ke + k ] += gbz[ (i*je + j)*ke + k ]*ez[ (i*je + j)*ke + k ];
            }   
        }   
    }
}

void calculate_fourier_transform_ex(int ie, int je, int ke, int number_of_frequencies, cvec& pt, rvec& ez, rvec& arg, double time_step, int kc){
    for (int j = 0; j < je; ++j){
        for (int i = 0; i < ie; ++i){
            for (int m = 0; m < number_of_frequencies; ++m){
                pt[ ((m*ie + i)*je + j)*ke + kc ] +=  complex<double>(cos(arg[m] * time_step), -sin(arg[m] * time_step)) * ez[(i*je + j)*ke + kc];
            }
        }
    }
}

void calculate_hx_field(int ie, int je, int ke, rvec& hx, rvec& ihx, rvec& ey, rvec& ez, rvec& fi1, rvec& fj2, rvec& fk2, rvec& fj3, rvec& fk3){
    #pragma omp for schedule(static) collapse(3)    
    for (int i = 0; i < ie; ++i){
        for (int j = 0; j < je - 1; ++j){
            for (int k = 0; k < ke - 1; ++k){
                double curl_e = ey[(i*je + j)*ke + (k + 1)] - ey[(i*je + j)*ke + k]  - ez[(i*je + (j + 1) )*ke + k] + ez[(i*je + j)*ke + k];
                ihx[(i*je + j)*ke + k] += curl_e;
                hx[(i*je + j)*ke + k] = fj3[j] * fk3[k] * hx[(i*je + j)*ke + k] + fj2[j] * fk2[k] * 0.5 * (curl_e + fi1[i] * ihx[(i*je + j)*ke + k] ); 
            }
        }
    }
}

void calculate_hy_field(int ie, int je, int ke, rvec& hy, rvec& ihy, rvec& ex, rvec& ez, rvec& fj1, rvec& fi2, rvec& fk2, rvec& fi3, rvec& fk3){
    #pragma omp for schedule(static) collapse(3)    
    for (int i = 0; i < ie - 1; ++i){
        for (int j = 0; j < je; ++j){
            for (int k = 0; k < ke - 1; ++k){
                double curl_e = ez[((i + 1)*je + j)*ke + k] - ez[( i*je + j)*ke + k]  - ex[(i*je + j )*ke + (k + 1) ] + ex[(i*je + j)*ke + k];
                ihy[(i*je + j)*ke + k] += curl_e;
                hy[(i*je + j)*ke + k] = fi3[i] * fk3[k] * hy[(i*je + j)*ke + k] + fi2[i] * fk2[k] * 0.5 * (curl_e + fj1[j] * ihy[(i*je + j)*ke + k] ); 
            }
        }
    }
}

void calculate_hz_field(int ie, int je, int ke, rvec& hz, rvec& ihz, rvec& ex, rvec& ey, rvec& fk1, rvec& fi2, rvec& fj2, rvec& fi3, rvec& fj3){
    #pragma omp for schedule(static) collapse(3)    
    for (int i = 0; i < ie - 1; ++i){
        for (int j = 0; j < je - 1; ++j){
            for (int k = 0; k < ke; ++k){
                double curl_e = ex[(i*je + (j + 1))*ke + k] - ex[( i*je + j)*ke + k]  - ey[((i + 1)*je + j )*ke + k ] + ey[(i*je + j)*ke + k];
                ihz[(i*je + j)*ke + k] += curl_e;
                hz[(i*je + j)*ke + k] = fi3[i] * fj3[j] * hz[(i*je + j)*ke + k] + fi2[i] * fj2[j] * 0.5 * (curl_e + fk1[k] * ihz[(i*je + j)*ke + k] ); 
            }
        }
    }
}

void calculate_hx_with_incident_field(int ie, int je, int ke, int ia, int ib, int ja, int jb, int ka, int kb, rvec& hx, rvec& ez_inc){
    #pragma omp for schedule(static) collapse(2)
    for (int i = ia; i <= ib; ++i){
        for (int k = ka; k <= kb; ++k){
            hx[( i*je + (ja - 1) )*ke + k] += 0.5*ez_inc[ja];
            hx[( i*je + jb)*ke + k] -= 0.5*ez_inc[jb]; 
        }    
    }
}


void calculate_hx_inc(int je, rvec& hx_inc, rvec& ez_inc){
    #pragma omp for schedule(static) 
    for (int j = 0; j < je - 1; ++j){
        hx_inc[j] += 0.5 * (ez_inc[j] - ez_inc[j + 1]);
    }
}

void calculate_hy_with_incident_field(int ie, int je, int ke, int ia, int ib, int ja, int jb, int ka, int kb, rvec& hy, rvec& ez_inc){
    #pragma omp for schedule(static) collapse(2)
    for (int j = ja; j <= jb; ++j){
        for (int k = ka; k <= kb; ++k){
            hy[( (ia - 1)*je + j)*ke + k] -= 0.5*ez_inc[j];
            hy[( ib*je + jb)*ke + k] += 0.5*ez_inc[j]; 
        }    
    }
}

void scalar_mult(rvec &input, rvec &output, const double scalar){
   
    size_t size = input.size();
    for (size_t i = 0; i < size; ++i){
        output[i] = input[i]*scalar;
    }
}

rvec calculate_field_strength(int ie, int je, int ke, rvec& fieldx, rvec& fieldy, rvec& fieldz){
    rvec sfield(ie*je*ke);

    #pragma omp for schedule(static) collapse(3)
    for (int i = 0; i < ie; ++i){
        for (int j = 0; j < je; ++j){
            for (int k = 0; k < ke; ++k){
                sfield[ (i*je + j)*ke + k ] = sqrt( fieldx[ (i*je + j)*ke + k ] * fieldx[(i*je + j)*ke + k] + fieldy[ (i*je + j)*ke + k ] * fieldy[(i*je + j)*ke + k] 
                                                    + fieldz[ (i*je + j)*ke + k ] * fieldz[(i*je + j)*ke + k]  );
            }
        }
    }
    return sfield;
}

void write_field(int ie, int je, int ke, rvec& fieldx, rvec& fieldy, rvec& fieldz, rvec& field_strength, string filename){
    ofstream field_file(filename);
    if (field_file.is_open()){
        field_file << "X(cm)" << "\t" << "Y(cm)" << "\t" << "Z (cm)" << "\t" <<"Field_X(x,y,z)" << "\t" <<  "Field_Y(x,y,z)" << "\t" << "Field_Z(x,y,z)" << "\t" << "Field_Strength(x,y,z)" << endl;
        for (int k = 0; k < ke; ++k){
            for (int j = 0; j < je; ++j){
                for (int i = 0; i < ie; ++i){
                    field_file << i << "\t" << j << " \t" << k << "\t" << fieldx[ (i*je + j)*ke + k ] << "\t" << fieldy[ (i*je + j)*ke + k ] << "\t" << fieldz[ (i*je + j)*ke + k ] << "\t" << field_strength[ (i*je + j)*ke + k ] << endl;
                }
            }
        }
    }
    else{
        cout << "Unable to open file" << endl;
    }
}

void write_amp_file(int je, rvec& amp, string filename){
    ofstream amp_file(filename);
    if (amp_file.is_open()){
        amp_file << "X (cm)" << "\t" << "Amplitud 50 MHz" << "\t" << "Amplitud 100 MHz" << "\t" << "Amplitud 200 MHz" << endl;
        for (int j = 0; j < je; ++j){
            amp_file << j - 20  << "\t" <<amp[j] << "\t" << amp[je + j] << "\t" << amp[2*je + j] << endl;
        }
    }
    else{
        cout << "Unable to open file" << endl;
    }
}


int main(){
    int ie, je, ke, ic, jc, kc, ia, ja, ka, ib, jb, kb;
    ie = 40;
    je = 40;
    ke = 40;
    ic = ie / 2;
    jc = je / 2;
    kc = ke / 2;
    ia = 7;
    ja = 7;
    ka = 7;
    ib = ie - ia -1;
    jb = je - ja -1;
    kb = ke - ka - 1;

    rvec ex(ie*je*ke, 0.0);
    rvec ey(ie*je*ke, 0.0);
    rvec ez(ie*je*ke, 0.0);
    rvec ix(ie*je*ke, 0.0);
    rvec iy(ie*je*ke, 0.0);
    rvec iz(ie*je*ke, 0.0);
    rvec dx(ie*je*ke, 0.0);
    rvec dy(ie*je*ke, 0.0);
    rvec dz(ie*je*ke, 0.0);
    rvec idx(ie*je*ke, 0.0);
    rvec idy(ie*je*ke, 0.0);
    rvec idz(ie*je*ke, 0.0);
    rvec hx(ie*je*ke, 0.0);
    rvec hy(ie*je*ke, 0.0);
    rvec hz(ie*je*ke, 0.0);
    rvec ihx(ie*je*ke, 0.0);
    rvec ihy(ie*je*ke, 0.0);
    rvec ihz(ie*je*ke, 0.0);
    
    rvec gax(ie*je*ke, 1.0);
    rvec gay(ie*je*ke, 1.0);
    rvec gaz(ie*je*ke, 1.0);
    rvec gbx(ie*je*ke, 0.0);
    rvec gby(ie*je*ke, 0.0);
    rvec gbz(ie*je*ke, 0.0);
    rvec hx_inc(je, 0.0);
    rvec ez_inc(je, 0.0);
    
    double ddx = 0.01;
    double dt = ddx / 6.0e8;
    double epsz = 8.854e-12;


    int number_of_frequencies = 3;
    rvec freq = {50.0e6, 200.0e6, 500.0e6};
    rvec arg(freq.size()); 
    scalar_mult(freq, arg, 2.0*PI*dt);
    cvec in(number_of_frequencies, complex<double>(0.0, 0.0) );
    cvec pt(number_of_frequencies*ie*je*ke, complex<double>(0.0, 0.0) );
    rvec amp(number_of_frequencies*je, 0.0);

    //specify dielectric sphere
    rvec epsilon(2, 1.0);
    rvec sigma(2, 0.0);
    epsilon[1] = 30.0;
    sigma[1] = 0.3;
    double radius = 10.0;
   
    #pragma omp for schedule(static) collapse(3)   
    for (int i = ia; i <= ib; ++i){
        for (int j = ja; j <= jb; ++j){
            for (int k = ka; k <= kb; ++k){
                double eps = epsilon[0];
                double cond = sigma[0];
                double xdist = static_cast<double>(ic) - static_cast<double>(i) - 0.5;
                double ydist = static_cast<double>(jc) - static_cast<double>(j);
                double zdist = static_cast<double>(kc) - static_cast<double>(k);
                double dist = sqrt( xdist * xdist + ydist * ydist + zdist * zdist);
                if (dist <= radius){
                    eps = epsilon[1];
                    cond = sigma[1]; 
                }
                gax[ (i*je + j)*ke + k ] = 1.0 / (eps + (cond * dt / epsz) );
                gbx[ (i*je + j)*ke + k ] = cond * dt / epsz;
            }
        }
    }

    #pragma omp for schedule(static) collapse(3)   
    for (int i = ia; i <= ib; ++i){
        for (int j = ja; j <= jb; ++j){
            for (int k = ka; k <= kb; ++k){
                double eps = epsilon[0];
                double cond = sigma[0];
                double xdist = static_cast<double>(ic) - static_cast<double>(i);
                double ydist = static_cast<double>(jc) - static_cast<double>(j) - 0.5;
                double zdist = static_cast<double>(kc) - static_cast<double>(k);
                double dist = sqrt( xdist * xdist + ydist * ydist + zdist * zdist);
                if (dist <= radius){
                    eps = epsilon[1];
                    cond = sigma[1];
                }
                gay[ (i*je + j)*ke + k ] = 1.0 / (eps + (cond * dt / epsz) );
                gby[ (i*je + j)*ke + k ] = cond * dt / epsz;
            }
        }
    }

    #pragma omp for schedule(static) collapse(3)   
    for (int i = ia; i <= ib; ++i){
        for (int j = ja; j <= jb; ++j){
            for (int k = ka; k <= kb; ++k){
                double eps = epsilon[0];
                double cond = sigma[0];
                double xdist = static_cast<double>(ic) - static_cast<double>(i);
                double ydist = static_cast<double>(jc) - static_cast<double>(j);
                double zdist = static_cast<double>(kc) - static_cast<double>(k) - 0.5;
                double dist = sqrt( xdist * xdist + ydist * ydist + zdist * zdist);
                if (dist <= radius){
                    eps = epsilon[1];
                    cond = sigma[1]; 
                }
                gaz[ (i*je + j)*ke + k ] = 1.0 / (eps + (cond * dt / epsz) );
                gbz[ (i*je + j)*ke + k ] = cond * dt / epsz;
            }
        }
    }

    //Pulse parameters
    double t0 = 20.0;
    double spread = 8.0;

    //Calculate the PML parameters
    int npml = 8;
    rvec gi1(ie, 0.0);
    rvec gi2(ie, 1.0);
    rvec gi3(ie, 1.0);
    rvec fi1(ie, 0.0);
    rvec fi2(ie, 1.0);
    rvec fi3(ie, 1.0);

    rvec gj1(je, 0.0);
    rvec gj2(je, 1.0);
    rvec gj3(je, 1.0);
    rvec fj1(je, 0.0);
    rvec fj2(je, 1.0);
    rvec fj3(je, 1.0);

    rvec gk1(ke, 0.0);
    rvec gk2(ke, 1.0);
    rvec gk3(ke, 1.0);
    rvec fk1(ke, 0.0);
    rvec fk2(ke, 1.0);
    rvec fk3(ke, 1.0);

    calculate_pml_parameters(npml, ie, je, ke, gi1, gi2, gi3, fi1, fi2, fi3, gj1, gj2, gj3, fj1, fj2, fj3, gk1, gk2, gk3, fk1, fk2, fk3);
    
    rvec boundary_low(2, 0.0);
    rvec boundary_high(2, 0.0);

    int nsteps = 500;

    //Main FDTD Loop
    for (int time_step = 1; time_step <= nsteps; ++time_step){
        for (int j = 1; j < je - 1; ++j){
            ez_inc[j] += 0.5* (hx_inc[j -1] - hx_inc[j]);
        }

        for (int m = 0; m < number_of_frequencies; ++m){
            in[m] += complex<double>(cos(arg[m] * static_cast<double>(time_step)) , -sin(arg[m] * static_cast<double>(time_step)) ) * ez_inc[ja - 1]; 
        }

        //AbsorbingBoundary Conditions
        double temp_low = boundary_low[0];
        boundary_low[0] = boundary_low[1];
        boundary_low[1] = ez_inc[1];
        ez_inc[0] = temp_low;

        double temp_high = boundary_high[0];
        boundary_high[0] = boundary_high[1];
        boundary_high[1] = ez_inc[je - 2];
        ez_inc[je - 1] = temp_high;
        
        //Calculate the D Fields
        calculate_dx_fied(ie, je, ke, dx, idx, hy, hz, gj3, gk3, gj2, gk2, gi1);
        calculate_dy_fied(ie, je, ke, dy, idy, hx, hz, gi3, gk3, gi2, gk2, gj1);
        calculate_dz_fied(ie, je, ke, dz, idz, hx, hy, gi3, gj3, gi2, gj2, gk1);
        
        //Add the source at the gap
        double pulse = exp(-0.5 * pow( (t0 - static_cast<double>(time_step)) / spread, 2) );
        ez_inc[3] = pulse;

        calculate_inc_dy_field(ie, je, ke, ia, ib, ja, jb, ka, kb, dy, hx_inc);
        calculate_inc_dz_field(ie, je, ke, ia, ib, ja, jb, ka, kb, dz, hx_inc);

        // Calculate the E field from the D field
        calculate_e_fields(ie, je, ke, dx, dy, dz, gax, gay, gaz, gbx, gby, gbz, ex, ey, ez, ix, iy, iz);
        
        //Calculate the Fourier transform of Ex
        calculate_fourier_transform_ex(ie, je, ke, number_of_frequencies, pt, ez, arg, static_cast<double>(time_step), kc);
        
        //Calculate the H fields
        calculate_hx_inc(je, hx_inc, ez_inc);
        calculate_hx_field(ie, je, ke, hx, ihx, ey, ez, fi1, fj2, fk2, fj3, fk3);
        calculate_hx_with_incident_field(ie, je, ke, ia, ib, ja, jb, ka, kb, hx, ez_inc);
        calculate_hy_field(ie, je, ke, hy, ihy, ex, ez, fj1, fi2, fk2, fi3, fk3);
        calculate_hy_with_incident_field(ie, je, ke, ia, ib, ja, jb, ka, kb, hy, ez_inc);
        calculate_hz_field(ie, je, ke, hz, ihz, ex, ey, fk1, fi2, fj2, fi3, fj3);

    }
    //Calculate the strenght o the electric_field
    rvec e_strength = calculate_field_strength(ie, je, ke, ex, ey, ez);
    rvec h_strength = calculate_field_strength(ie, je, ke, hx, hy, hz);

    rvec amp_in(number_of_frequencies);
    for (int m = 0; m < number_of_frequencies; ++m) {
        amp_in[m] = abs(in[m]);
    }

    // Calculate the Fourier amplituyde of the total field
    for (int m = 0; m < number_of_frequencies; ++m){
        for (int j = ja; j <= jb; ++j){
            if (gaz[ (ic*je +j)*ke + kc ] < 1){
                amp[ m*je + j ] = 1 / (amp_in[m]) * abs(pt[ ((m*ie + ic)*je + j)*ke + kc ]); 
            }
        }
    }

    write_amp_file(je, amp, "../data/amp_data.txt");
    write_field(ie, je, ke, ex, ey, ez, e_strength, "../data/E_field.txt");
    write_field(ie, je, ke, hx, hy, hz, h_strength, "../data/H_field.txt");
    
    return 0;
}
