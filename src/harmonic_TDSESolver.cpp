#include "propag/config.h"
#include "propag/harmonic_TDSESolver.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

void createDirectories(std::string name){
    try{
        bool status = std::filesystem::create_directories(name);
        if (!status)
            throw(status);
    }
    catch(bool status){
        SPDLOG_ERROR("[HarmonicTDSESolver] Temporal directories could not be created.");
    }
}

void writeField(Field *field, int n, std::string tmppath){
    std::string path_efieldx = tmppath+"/efieldx.dat";
    std::ofstream efieldx(path_efieldx);
    if (efieldx.is_open()){
        for(int i = 0; i<n; i++){
            efieldx << std::setprecision(12) << std::fixed << field->mFieldX[i].real()<<std::endl;
        }
        efieldx.close();
    }

    std::string path_efieldy = tmppath+"/efieldy.dat";
    std::ofstream efieldy(path_efieldy);
    if (efieldy.is_open()){
        for(int i = 0; i<n; i++){
            efieldy << std::setprecision(12) << std::fixed << field->mFieldY[i].real()<<std::endl;
        }
        efieldy.close();
    }

    std::string path_efieldz = tmppath+"/efieldz.dat";
    std::ofstream efieldz(path_efieldz);
    if (efieldz.is_open()){
        for(int i = 0; i<n; i++){
            efieldz << std::setprecision(12) << std::fixed << field->mFieldZ[i].real()<<std::endl;
        }
        efieldz.close();
    }
}

void readAccel(std::string tmppath, double *accel_x, double *accel_y, double *accel_z, const int n){
    std::string path_accel_x = tmppath + "/acc_i.dat";
    std::string path_accel_y = tmppath + "/acc_j.dat";
    std::string path_accel_z = tmppath + "/acc_k.dat";
    
    std::ifstream file_accel_x(path_accel_x);
    std::complex<double> value = (0.0,0.0);
    if(file_accel_x.is_open()){
        for(int i=0; i<n; i++){
            file_accel_x >> value;
            accel_x[i] = value.real();
        }
        file_accel_x.close();
    }

    std::ifstream file_accel_y(path_accel_y);
    if(file_accel_y.is_open()){
        for(int i=0; i<n; i++){
            file_accel_y >> value;
            accel_y[i] = value.real();
        }
        file_accel_y.close();
    }
/*
    std::ifstream file_accel_z(path_accel_z);
    if(file_accel_z.is_open()){
        for(int i=0; i<n; i++){
            file_accel_z >> value;
            accel_z[i] = value.real();
        }
        file_accel_z.close();
    }
*/
}

HarmonicTDSESolver::HarmonicTDSESolver(std::shared_ptr<Settings> settings) : Harmonic(settings){
    param = new Parameters();

    param->n_threads       =  mSettings->mData.numThreads;
    param->use_potential  =  1;
    param->init_wf        =  mSettings->mData.tdsesolver_initwf_type;
    param->geometry       =  mSettings->mData.tdsesolver_coordinateType;
    param->ni             =  mSettings->mData.tdsesolver_ndx;
    param->nj             =  mSettings->mData.tdsesolver_ndy;
    param->nk             =  mSettings->mData.tdsesolver_ndz;
    param->imin           = -mSettings->mData.tdsesolver_xmax;
    param->imax           =  mSettings->mData.tdsesolver_xmax;
    param->jmin           = -mSettings->mData.tdsesolver_ymax;
    param->jmax           =  mSettings->mData.tdsesolver_ymax;
    param->kmin           = -mSettings->mData.tdsesolver_zmax;
    param->kmax           =  mSettings->mData.tdsesolver_zmax;
    param->tmax_ev        =  mSettings->TMAX;
    param->tmax_sim       =  mSettings->TMAX;
    param->nt             =  mSettings->NDT;
    param->dt             =  mSettings->DT;
    param->dt_ITP         =  mSettings->mData.tdsesolver_dtITP;
    param->nt_ITP         =  mSettings->mData.tdsesolver_ntITP;
    param->nt_diag        =  mSettings->mData.tdsesolver_ntDiag;
        
    param->check_param();
    param->print();
}

HarmonicTDSESolver::~HarmonicTDSESolver(){
    delete param;
    delete tdsesolver;
}

void HarmonicTDSESolver::calculateAcceleration(Field *field){
    std::string timestamp = std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                          std::chrono::system_clock::now().time_since_epoch()).count());
    std::string tmp_path_results = "tmp/results_" + timestamp;
    std::string tmp_path_efields = "tmp/efields_" + timestamp;
    createDirectories(tmp_path_results);
    createDirectories(tmp_path_efields);
    param->use_field_file =  1;
    param->file_fieldx    =  tmp_path_efields + "/efieldx.dat";
    param->file_fieldy    =  tmp_path_efields + "/efieldy.dat";
    param->file_fieldz    =  tmp_path_efields + "/efieldz.dat";
    param->path_results   =  tmp_path_results;
    param->n_probes       =  3;
    param->probe_def      =  "acc_i," + tmp_path_results + "/acc_i.dat;"
                             "acc_j," + tmp_path_results + "/acc_j.dat;"
                             "acc_k," + tmp_path_results + "/acc_k.dat";


    writeField(field, param->nt, tmp_path_efields);
    double *accel_z = nullptr;
    tdsesolver = new TDSESolver(param);
    tdsesolver->ipropagate();
    tdsesolver->propagate();
    readAccel(tmp_path_results, mAccelX, mAccelY, accel_z, param->nt);
    delete tdsesolver;
}
