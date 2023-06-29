#include "data_loader.h"
#include "hyperef.h"
bool verify(std::string file_l,std::string file_r);
int main(int argc, char **argv) {
	
	// parse command-line arguments
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " <dataset>" << std::endl;
    return 0;
  }
  std::string filename = argv[1];
  printf("%zu  %zu\n", sizeof(float),sizeof(float_t));
  printf("%zu  %zu\n", sizeof(double),sizeof(double_t));
  // std::string filename = "data/ibm05.hgr";
  spmv::io::CSCMatrix<val_t> ar;
  ar = load_csr_matrix_from_float_hmetis_unweighted(filename);
  uint32_t L, R;
  L = 1;
  R = 1;
  HyperEF(ar,L,R);
/*
  if(!verify("Debug","../hyperef/Debug")){
      std::cout << "FATAL:: Filter C++ implementation do not match!" <<std::endl;
  }*/

  return 0;
}

bool verify(std::string file_l,std::string file_r) {
    bool status;
    float epsilon = 0.0001;
    std::vector<double> rv;
    std::string line;
    std::ifstream myFile;            //creates stream myFile
    myFile.open(file_l.c_str());  //opens .txt file
    uint32_t lines_read1;
    uint32_t lines_read2;
    double residue;
    lines_read1=0;
    lines_read2=0;
    residue=0.0;
    status=true;
    while (myFile.good()) {
        line.clear();
        getline(myFile, line);
        //rv[lines_read1++]=std::strtof(line.c_str(),nullptr);
        rv.push_back(std::strtof(line.c_str(),nullptr));
        lines_read1++;
    }
    myFile.close();
    myFile.open(file_r.c_str());
    while (myFile.good()) {
        line.clear();
        getline(myFile, line);
        if(lines_read2==lines_read1){
            std::cout << "Length mismatch" << std::endl;
            status=false;
            break;
        }
        rv[lines_read2]-=std::strtof(line.c_str(),nullptr);
        residue+=abs(rv[lines_read2]);
        if(rv[lines_read2++]>epsilon){
            std::cout << "Error: Result mismatch"
                      << std::endl;
            status=false;
            break;
        }
    }
    if(residue > pow(epsilon,0.75))
        status=false;
    myFile.close();
    return status;
}