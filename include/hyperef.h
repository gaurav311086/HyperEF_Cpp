#include <math.h>
#include <random>
#include <algorithm>
#include <numeric>
#include "data_loader.h"

typedef  std::double_t val_t;

// Data structure for csc matrix.
template<typename data_type>
struct SVEC {
    /*! \brief The number of rows of the sparse vector */
    uint32_t size;
    /*! \brief The non-zero data of the sparse vector */
    std::vector<data_type> adj_data;
    /*! \brief The index pointers of the sparse vector */
    std::vector<uint32_t> adj_indptr;
};

/*int get_random(int min, int max)
{
    static std::mt19937 mt{ std::random_device{}() };
    std::uniform_int_distribution die{ min, max };

    return die(mt);
}*/
template<typename data_type>
void dump_var(std::vector <data_type> var, std::string file){
    return;
    std::ofstream oFile;
    oFile.open(file.c_str());
    oFile.precision(19);
    for(int i =0; i < var.size();i++){
        oFile << var[i] << std::endl;
    }
    oFile.close();
}


// Load a csc matrix from a hmetis (unweighted hypergraph) file. The sparse matrix should have float data type.
spmv::io::CSCMatrix<val_t> load_csr_matrix_from_float_hmetis_unweighted(std::string csc_float_hmetis_path) {
    spmv::io::CSCMatrix<val_t> csc_matrix;
    std::ifstream io;
    io.open(csc_float_hmetis_path.c_str());
    if ( !io.is_open() ){
        throw std::runtime_error("FATAL: failed ifstream open");
    }
    std::string line;
    uint32_t lines_read;
    uint32_t this_hedge_node;
    uint32_t nnz_elements_read;
    nnz_elements_read = 0;
    lines_read=0;
    while (io.good()){
        bool first_int;
        bool second_int;
        line.clear();
        getline(io, line);
        std::string substr;
        substr.clear();
        first_int = false;
        second_int = false;
        for(uint32_t i=0; i < line.length();i++){
          if(line[i] == ' ' || (i==line.length()-1) ) {
              if(!(line[i]=='\r' || line[i]=='\n' || line[i]==' ')) {
                  substr.push_back(line[i]);
              }
              if(substr.length() > 0){
                if(lines_read>0){
                    nnz_elements_read++;
                    this_hedge_node = std::stoi(substr.c_str());
                  csc_matrix.adj_indices.push_back(this_hedge_node-1);
                  csc_matrix.adj_data.push_back(1.0);
                } else {
                  if(first_int){
                      if(second_int){
                          io.close();
                          throw std::runtime_error("FATAL: Only unweighted hmetis files supported");
                      }
                      csc_matrix.num_rows = std::stoi(substr.c_str());
                      second_int=true;
                  } else {
                    csc_matrix.num_cols = std::stoi(substr.c_str());
                    first_int = true;
                  }
                }
              }
              substr.clear();
          } else {
              if(!(line[i]=='\r' || line[i]=='\n')) {
                      substr.push_back(line[i]);
              } else {
                substr.clear();
              }
          }
        }
        if(line.length()>0) {
            csc_matrix.adj_indptr.push_back(nnz_elements_read);
        }
        lines_read++;
    }
    io.close();
    return csc_matrix;
}

// Convert csc to csr.
template<typename data_type>
spmv::io::CSRMatrix<data_type> csc2csr(spmv::io::CSCMatrix<data_type> const &csc_matrix) {
    spmv::io::CSRMatrix<data_type> csr_matrix;
    csr_matrix.num_rows = csc_matrix.num_rows;
    csr_matrix.num_cols = csc_matrix.num_cols;
    csr_matrix.adj_data = std::vector<data_type>(csc_matrix.adj_data.size());
    csr_matrix.adj_indices = std::vector<uint32_t>(csc_matrix.adj_indices.size());
    csr_matrix.adj_indptr = std::vector<uint32_t>(csr_matrix.num_rows + 1);
    // Convert adj_indptr
    uint32_t nnz = csc_matrix.adj_indptr[csc_matrix.num_cols];
    std::vector<uint32_t> nnz_each_row(csr_matrix.num_rows);
    std::fill(nnz_each_row.begin(), nnz_each_row.end(), 0);
    for (size_t n = 0; n < nnz; n++) {
        nnz_each_row[csc_matrix.adj_indices[n]]++;
    }
    csr_matrix.adj_indptr[0] = 0;
    for (size_t row_idx = 0; row_idx < csr_matrix.num_rows; row_idx++) {
        csr_matrix.adj_indptr[row_idx + 1] = csr_matrix.adj_indptr[row_idx] + nnz_each_row[row_idx];
    }
    assert(csr_matrix.adj_indptr[csr_matrix.num_rows] == nnz);
    // Convert adj_data and adj_indices
    std::vector<uint32_t> nnz_consumed_each_row(csr_matrix.num_rows);
    std::fill(nnz_consumed_each_row.begin(), nnz_consumed_each_row.end(), 0);
    for (size_t col_idx = 0; col_idx < csc_matrix.num_cols; col_idx++){
        for (size_t i = csc_matrix.adj_indptr[col_idx]; i < csc_matrix.adj_indptr[col_idx + 1]; i++){
            uint32_t row_idx = csc_matrix.adj_indices[i];
            uint32_t dest = csr_matrix.adj_indptr[row_idx] + nnz_consumed_each_row[row_idx];
            csr_matrix.adj_indices[dest] = col_idx;
            csr_matrix.adj_data[dest] = csc_matrix.adj_data[i];
            nnz_consumed_each_row[row_idx]++;
        }
    }
    for (size_t row_idx = 0; row_idx < csr_matrix.num_rows; row_idx++) {
        assert(nnz_consumed_each_row[row_idx] == nnz_each_row[row_idx]);
    }
    return csr_matrix;
}

template<typename data_type>
spmv::io::CSCMatrix<data_type> transpose(spmv::io::CSCMatrix<data_type> const &csc_matrix) {
    spmv::io::CSCMatrix<data_type> csr_matrix;
    csr_matrix.num_rows = csc_matrix.num_cols;
    csr_matrix.num_cols = csc_matrix.num_rows;
    csr_matrix.adj_data = std::vector<data_type>(csc_matrix.adj_data.size());
    csr_matrix.adj_indices = std::vector<uint32_t>(csc_matrix.adj_indices.size());
    csr_matrix.adj_indptr = std::vector<uint32_t>(csr_matrix.num_rows + 1);
    // Convert adj_indptr
    uint32_t nnz = csc_matrix.adj_indptr[csc_matrix.num_cols];
    std::vector<uint32_t> nnz_each_row(csr_matrix.num_rows);
    std::fill(nnz_each_row.begin(), nnz_each_row.end(), 0);
    for (size_t n = 0; n < nnz; n++) {
        nnz_each_row[csc_matrix.adj_indices[n]]++;
    }
    csr_matrix.adj_indptr[0] = 0;
    for (size_t row_idx = 0; row_idx < csr_matrix.num_rows; row_idx++) {
        csr_matrix.adj_indptr[row_idx + 1] = csr_matrix.adj_indptr[row_idx] + nnz_each_row[row_idx];
    }
    assert(csr_matrix.adj_indptr[csr_matrix.num_rows] == nnz);
    // Convert adj_data and adj_indices
    std::vector<uint32_t> nnz_consumed_each_row(csr_matrix.num_rows);
    std::fill(nnz_consumed_each_row.begin(), nnz_consumed_each_row.end(), 0);
    for (size_t col_idx = 0; col_idx < csc_matrix.num_cols; col_idx++){
        for (size_t i = csc_matrix.adj_indptr[col_idx]; i < csc_matrix.adj_indptr[col_idx + 1]; i++){
            uint32_t row_idx = csc_matrix.adj_indices[i];
            uint32_t dest = csr_matrix.adj_indptr[row_idx] + nnz_consumed_each_row[row_idx];
            csr_matrix.adj_indices[dest] = col_idx;
            csr_matrix.adj_data[dest] = csc_matrix.adj_data[i];
            nnz_consumed_each_row[row_idx]++;
        }
    }
    for (size_t row_idx = 0; row_idx < csr_matrix.num_rows; row_idx++) {
        assert(nnz_consumed_each_row[row_idx] == nnz_each_row[row_idx]);
    }
    return csr_matrix;
}

template<typename data_type>
spmv::io::CSCMatrix<data_type> create_diagonal_matrix(std::vector<data_type> const &W, uint32_t rows_cols = 0) {
    spmv::io::CSCMatrix<data_type> csc_matrix;
    uint32_t diagonal_elements;
    uint32_t nWeights;
    if(rows_cols)
        diagonal_elements = rows_cols;
    else
        diagonal_elements = W.size();
    nWeights = W.size();
    csc_matrix.num_rows=diagonal_elements;
    csc_matrix.num_cols=diagonal_elements;
    csc_matrix.adj_data = std::vector<data_type>(diagonal_elements);
    csc_matrix.adj_indices = std::vector<uint32_t>(diagonal_elements);
    csc_matrix.adj_indptr = std::vector<uint32_t>(diagonal_elements+1);
    for(uint32_t i=0; i<diagonal_elements;i++){
        if(nWeights==1)
            csc_matrix.adj_data[i]=W[0];
        else
            csc_matrix.adj_data[i]=W[i];
        csc_matrix.adj_indices[i]=i;
        csc_matrix.adj_indptr[i]=i;
    }
    csc_matrix.adj_indptr[diagonal_elements]=diagonal_elements;
    return csc_matrix;
}

template<typename data_type>
spmv::io::CSCMatrix<data_type> sum(spmv::io::CSCMatrix<data_type> const &csc_matrix1, spmv::io::CSCMatrix<data_type> const &csc_matrix2){
    spmv::io::CSCMatrix<data_type> csc_matrixOut;
    uint32_t sz1;
    uint32_t sz2;
    uint32_t elemCovered1;
    uint32_t elemCovered2;
    uint32_t elemAdded;
    if(csc_matrix1.num_rows != csc_matrix2.num_rows)
        throw std::runtime_error("FATAL: Matrix dimensions do match..check rows");
    if(csc_matrix1.num_cols != csc_matrix2.num_cols)
        throw std::runtime_error("FATAL: Matrix dimensions do match..check columns");
    sz1 = csc_matrix1.adj_data.size();
    sz2 = csc_matrix2.adj_data.size();
    csc_matrixOut.num_rows=csc_matrix1.num_rows;
    csc_matrixOut.num_cols=csc_matrix2.num_cols;
    csc_matrixOut.adj_indptr = std::vector<uint32_t>(csc_matrixOut.num_cols + 1);
    elemCovered1 = 0;
    elemCovered2 = 0;
    elemAdded    = 0;

    csc_matrixOut.adj_indptr[0]=elemAdded;
    for(uint32_t col=0;col<csc_matrixOut.num_cols;col++){
        uint32_t nnz_ele_in_mat1;
        uint32_t nnz_ele_in_mat2;
        uint32_t ToBeTraversed;
        uint32_t nnz_traversed_this_step;
        nnz_ele_in_mat1 = csc_matrix1.adj_indptr[col+1]-csc_matrix1.adj_indptr[col];
        nnz_ele_in_mat2 = csc_matrix2.adj_indptr[col+1]-csc_matrix2.adj_indptr[col];
        nnz_traversed_this_step=0;
        ToBeTraversed = nnz_ele_in_mat1 + nnz_ele_in_mat2;
        for(uint32_t traversed=0; traversed < ToBeTraversed; traversed+=nnz_traversed_this_step){
            nnz_traversed_this_step=0;
            if((nnz_ele_in_mat1 == 0) && (nnz_ele_in_mat2==0)){
                nnz_traversed_this_step=0;
                //throw std::runtime_error("FATAL: Loop error");
            } else {
                if(nnz_ele_in_mat1 == 0) {
                    //matrix1 has no element in this row
                    //  but mat2 has some. copy those
                    uint32_t copy_ops;
                    copy_ops= nnz_ele_in_mat2;
                    for(uint32_t iter=0;iter<copy_ops;iter++){
                        csc_matrixOut.adj_data.push_back(csc_matrix2.adj_data[elemCovered2]);
                        csc_matrixOut.adj_indices.push_back(csc_matrix2.adj_indices[elemCovered2]);
                        elemCovered2++;
                        nnz_traversed_this_step++;
                        elemAdded++;
                        nnz_ele_in_mat2--;
                    }
                } else if(nnz_ele_in_mat2 == 0) {
                    //matrix2 has no element in this row
                    //  but mat1 has some. copy those
                    uint32_t copy_ops;
                    copy_ops= nnz_ele_in_mat1;
                    for(uint32_t iter=0;iter<copy_ops;iter++){
                        csc_matrixOut.adj_data.push_back(csc_matrix1.adj_data[elemCovered1]);
                        csc_matrixOut.adj_indices.push_back(csc_matrix1.adj_indices[elemCovered1]);
                        elemCovered1++;
                        nnz_traversed_this_step++;
                        elemAdded++;
                        nnz_ele_in_mat1--;
                    }
                } else {
                    //Both Matrix have some elements.
                    if(csc_matrix1.adj_indices[elemCovered1] == csc_matrix2.adj_indices[elemCovered2]){
                        // same row entry present..add
                        csc_matrixOut.adj_data.push_back(csc_matrix1.adj_data[elemCovered1]+csc_matrix2.adj_data[elemCovered2]);
                        csc_matrixOut.adj_indices.push_back(csc_matrix1.adj_indices[elemCovered1]);
                        elemCovered1++;
                        elemCovered2++;
                        nnz_traversed_this_step+=2;
                        nnz_ele_in_mat1--;
                        nnz_ele_in_mat2--;
                    } else if(csc_matrix1.adj_indices[elemCovered1] < csc_matrix2.adj_indices[elemCovered2]){
                        // append mat1 entry now
                        csc_matrixOut.adj_data.push_back(csc_matrix1.adj_data[elemCovered1]);
                        csc_matrixOut.adj_indices.push_back(csc_matrix1.adj_indices[elemCovered1]);
                        elemCovered1++;
                        nnz_traversed_this_step+=1;
                        nnz_ele_in_mat1--;
                    } else {
                        //append mat2 entry now
                        csc_matrixOut.adj_data.push_back(csc_matrix2.adj_data[elemCovered2]);
                        csc_matrixOut.adj_indices.push_back(csc_matrix2.adj_indices[elemCovered2]);
                        elemCovered2++;
                        nnz_traversed_this_step+=1;
                        nnz_ele_in_mat2--;
                    }
                    elemAdded++;
                }
            }
        }
        csc_matrixOut.adj_indptr[col+1]=elemAdded;
    }
    return csc_matrixOut;
}

spmv::io::CSCMatrix<val_t> StarW(spmv::io::CSCMatrix<val_t> const &ar, std::vector<val_t> &W){
    uint32_t mx;
    uint32_t sz;

    mx = ar.num_rows;
    sz = ar.num_cols;

    spmv::io::CSCMatrix <val_t> matrix;
    spmv::io::CSCMatrix <val_t> matrixT;
    spmv::io::CSCMatrix <val_t> matrixO;
    matrix.num_rows= mx + sz;
    matrix.num_cols= mx + sz;
    matrixT.num_rows= mx + sz;
    matrixT.num_cols= mx + sz;
    matrix.adj_data = std::vector<val_t>(ar.adj_data.size());
    matrix.adj_indices = std::vector<uint32_t>(ar.adj_indices.size());
    matrix.adj_indptr = std::vector<uint32_t>(mx+sz + 1);
    for(uint32_t i=0; i<mx+1;i++){
        matrix.adj_indptr[i] = 0;
    }
    uint32_t indices_covered;
    uint32_t lastArIndptr;
    indices_covered = 0;
    lastArIndptr = ar.adj_indptr[0];
    for(uint32_t iter =0;iter<sz;iter++){
        uint32_t LN;
        val_t edge_weight;
        uint32_t currArIndptr;
        currArIndptr = ar.adj_indptr[iter+1];
        matrix.adj_indptr[iter+mx+1] = currArIndptr;
        LN= currArIndptr - lastArIndptr;
        edge_weight=1.0/val_t(LN);
        for(uint32_t j=0; j< LN; j++){
            matrix.adj_data[indices_covered]=edge_weight;
            matrix.adj_indices[indices_covered]=ar.adj_indices[indices_covered];
            indices_covered++;
        }
        lastArIndptr=currArIndptr;
    }
    matrixT = transpose(matrix);
    matrixO = sum(matrix,matrixT);
    return matrixO;

}

spmv::io::CSCMatrix<val_t> SubRoutine001Filter(spmv::io::CSCMatrix <val_t> &AD){
    // = sparse(I2, I2, sparsevec(dg))
    spmv::io::CSCMatrix<val_t> csc_matrix;
    uint32_t ElementsSummed;
    uint32_t ElementsAdded;
    csc_matrix.num_rows=AD.num_rows;
    csc_matrix.num_cols=AD.num_cols;
    csc_matrix.adj_indptr=std::vector<uint32_t>(csc_matrix.num_cols+1);
    ElementsSummed=AD.adj_indptr[0];
    ElementsAdded=0;
    for(uint32_t col=0;col<AD.num_cols;col++){
        uint32_t ElementsToBeSummed;
        val_t sum;
        ElementsToBeSummed = AD.adj_indptr[col+1] - ElementsSummed;
        sum=0;
        for(uint32_t ele=0;ele<ElementsToBeSummed;ele++){
            sum+=AD.adj_data[ele+ElementsSummed];
        }
        ElementsSummed = AD.adj_indptr[col+1];
        if(ElementsToBeSummed){
            ElementsAdded++;
            sum = pow(sum, (-0.5));
            csc_matrix.adj_data.push_back(sum);
            csc_matrix.adj_indices.push_back(col);
        }
        csc_matrix.adj_indptr[col+1]=ElementsAdded;
    }
    return csc_matrix;
}
std::vector<val_t> SubRoutine002Filter(std::vector <val_t> &rv){
    std::vector<val_t> new_vector;
    //=rv - ((dot(rv, on) / dot(on, on)) * on)
    uint32_t sz;
    val_t dot_product;
    val_t mean_dot_product;

    sz=rv.size();
    new_vector=std::vector<val_t>(sz);
    dot_product=0;
    for(uint32_t i=0;i<sz;i++) {
        dot_product += rv[i];
    }
    mean_dot_product=dot_product/val_t(sz);
    for(uint32_t i=0;i<sz;i++){
        new_vector[i]=rv[i]-mean_dot_product;
    }
    return new_vector;
}



std::vector<val_t> SubRoutine003Filter(std::vector <val_t> &sm_ot,val_t power=2){
    std::vector<val_t> new_vector;
    //= sm_ot ./ norm(sm_ot)
    uint32_t sz;
    val_t dot_product;
    val_t norm_dot_product;

    sz=sm_ot.size();
    new_vector=std::vector<val_t>(sz);
    dot_product=0;
    for(uint32_t i=0;i<sz;i++) {
        dot_product += sm_ot[i]*sm_ot[i];
    }
    norm_dot_product=sqrt(dot_product);
    for(uint32_t i=0;i<sz;i++){
        new_vector[i]=sm_ot[i]/norm_dot_product;
    }
    return new_vector;
}

std::vector<val_t> SubRoutine004Filter(spmv::io::CSCMatrix <val_t> &D, std::vector <val_t> &sm){
    //=D*sm
    std::vector<val_t> new_vector;
    uint32_t sz;
    sz=D.num_rows;
    new_vector=std::vector<val_t>(sz);
    std::fill(new_vector.begin(), new_vector.end(), 0);
    for (uint32_t row_idx = 0; row_idx < sz; row_idx++) {
        uint32_t start;
        uint32_t end;
        start = D.adj_indptr[row_idx];
        end = D.adj_indptr[row_idx + 1];
        for (uint32_t i = start; i < end; i++) {
            uint32_t idx;
            idx = D.adj_indices[i];
            new_vector[row_idx] += D.adj_data[i] * sm[idx];
        }
    }
    return new_vector;
}

std::vector<std::vector<val_t>>Filter_old(std::vector<val_t> rv, uint32_t k, spmv::io::CSCMatrix <val_t> &AD, uint32_t mx, uint32_t initial, uint32_t interval, uint32_t Ntot){
    std::vector<std::vector<val_t>> V(mx,std::vector<val_t>(Ntot));
    spmv::io::CSCMatrix<val_t> AD_diagnal;
    spmv::io::CSCMatrix<val_t> D;
    std::vector<val_t> AD_diagnal_val = {0.1};
    uint32_t sz;
    std::vector<std::vector<val_t>> sm_vec( mx , std::vector<val_t> (k));
    std::vector<val_t> sm_ot;
    std::vector<val_t> sm;
    std::vector<val_t> sm_norm;
    uint32_t count;

    sz = AD.num_rows;
    AD_diagnal = create_diagonal_matrix(AD_diagnal_val,sz);


    //V = zeros(mx, Ntot);

    //sm_vec = zeros(mx, k);

    //AD = AD .* 1.0

    //AD[diagind(AD, 0)] = AD[diagind(AD, 0)] .+ 0.1
    AD = sum(AD,AD_diagnal);

    //dg = sum(AD, dims = 1) .^ (-.5)

    //I2 = 1:sz

    //D = sparse(I2, I2, sparsevec(dg))
    D=SubRoutine001Filter(AD);

    //on = ones(Int, length(rv))

    //sm_ot = rv - ((dot(rv, on) / dot(on, on)) * on)
    sm_ot = SubRoutine002Filter(rv);

    //sm = sm_ot ./ norm(sm_ot);
    sm = SubRoutine003Filter(sm_ot);

    count = 1;

    for(uint32_t loop=0;loop<k;loop++){
        //sm = D * sm
        sm = SubRoutine004Filter(D,sm);

        //sm = AD * sm
        sm = SubRoutine004Filter(AD,sm);

        //sm = D * sm
        sm = SubRoutine004Filter(D,sm);

        //sm_ot = sm - ((dot(sm, on) / dot(on, on)) * on)
        sm_ot = SubRoutine002Filter(sm);

        //sm_norm = sm_ot ./ norm(sm_ot);
        sm_norm = SubRoutine003Filter(sm_ot);
        for(uint32_t j=0;j<mx;j++){
            sm_vec[j][ loop] = sm_norm[j];
        }

    }

    //V = sm_vec[:, interval:interval:end]
    for(uint32_t j=0;j<mx;j++){
        for(uint32_t i=0;i<Ntot;i++){
            V[j][i]=sm_vec[j][(i*interval)+interval-1];
        }
    }

    return V;
}

std::vector<std::vector<val_t>>Filter(std::vector<val_t> rv, uint32_t k, spmv::io::CSCMatrix <val_t> &AD, uint32_t mx, uint32_t initial, uint32_t interval, uint32_t Ntot){
    std::vector<std::vector<val_t>> V(mx,std::vector<val_t>(Ntot));
    spmv::io::CSCMatrix<val_t> AD_diagnal;
    spmv::io::CSCMatrix<val_t> D;
    std::vector<val_t> AD_diagnal_val = {0.1};
    uint32_t sz;
    std::vector<std::vector<val_t>> sm_vec( mx , std::vector<val_t> (k));
    std::vector<std::vector<val_t>> sm( k+1 , std::vector<val_t> (rv.size()));
    std::vector<val_t> ping(rv.size());
    std::vector<val_t> pong(rv.size());
    std::vector<val_t> sm_ot;
    //std::vector<val_t> sm;
    std::vector<val_t> sm_norm;
    uint32_t count;

    sz = AD.num_rows;
    AD_diagnal = create_diagonal_matrix(AD_diagnal_val,sz);


    //V = zeros(mx, Ntot);

    //sm_vec = zeros(mx, k);

    //AD = AD .* 1.0

    //AD[diagind(AD, 0)] = AD[diagind(AD, 0)] .+ 0.1
    AD = sum(AD,AD_diagnal);

    //dg = sum(AD, dims = 1) .^ (-.5)

    //I2 = 1:sz

    //D = sparse(I2, I2, sparsevec(dg))
    D=SubRoutine001Filter(AD);

    //on = ones(Int, length(rv))

    //sm_ot = rv - ((dot(rv, on) / dot(on, on)) * on)
    sm_ot = SubRoutine002Filter(rv);

    //sm = sm_ot ./ norm(sm_ot);
    sm[0] = SubRoutine003Filter(sm_ot);

    count = 1;

    std::string dfile;
    for(uint32_t loop=0;loop<k;loop++){
        //sm = D * sm
        pong = SubRoutine004Filter(D,sm[loop]);
        dfile="sm_op1_"+std::to_string(loop);
        dump_var(sm[loop],dfile);

        //sm = AD * sm
        ping = SubRoutine004Filter(AD,pong);
        dfile="sm_op2_"+std::to_string(loop);
        dump_var(pong,dfile);

        //sm = D * sm
        sm[loop+1] = SubRoutine004Filter(D,ping);
        dfile="sm_op3_"+std::to_string(loop);
        dump_var(ping,dfile);

    }
    for(uint32_t loop=0;loop<k;loop++){

        //sm_ot = sm - ((dot(sm, on) / dot(on, on)) * on)
        sm_ot = SubRoutine002Filter(sm[loop+1]);
        dfile="sm_op4_"+std::to_string(loop);
        dump_var(sm_ot,dfile);
        //sm_norm = sm_ot ./ norm(sm_ot);
        sm_norm = SubRoutine003Filter(sm_ot);
        dfile="sm_op5_"+std::to_string(loop);
        dump_var(sm_norm,dfile);
        for(uint32_t j=0;j<mx;j++){
            sm_vec[j][ loop] = sm_norm[j];
        }
    }


    //V = sm_vec[:, interval:interval:end]
    for(uint32_t j=0;j<mx;j++){
        for(uint32_t i=0;i<Ntot;i++){
            V[j][i]=sm_vec[j][(i*interval)+interval-1];
        }
    }

    return V;
}


std::vector<std::vector<val_t>> QR_decomposition_return_Q_only(std::vector<std::vector<val_t>> const &SV){
    uint32_t  SV_rows;
    uint32_t  SV_cols;
    SV_rows = SV.size();
    SV_cols = SV[0].size();

    std::vector<std::vector<val_t>> Q_t(SV_cols, std::vector<val_t> (SV_rows));
    std::vector<std::vector<val_t>> q(SV_rows, std::vector<val_t> (SV_cols));
    for(uint32_t col=0;col<SV_cols;col++){
        for(uint32_t row=0;row<SV_rows;row++) {
            Q_t[col][row]=SV[row][col];
        }
    }
    val_t uu,up;
    for(uint32_t col=0;col<SV_cols;col++){
        //compute vector
        std::vector<val_t> aCol;
        aCol=Q_t[col];
        for(int32_t col_proj=0;col_proj<col;col_proj++) {
            up=0;
            for(uint32_t row=0;row<SV_cols;row++) {
                up+=aCol[row]*Q_t[col_proj][row];
            }
            //up=up/1
            for(uint32_t row=0;row<SV_cols;row++) {
                Q_t[col][row]-=Q_t[col_proj][row]*up;
            }
        }
        //normalize vector
        uu=0;
        for(uint32_t row=0;row<SV_cols;row++) {
            uu+=Q_t[col][row]*Q_t[col][row];
        }
        uu=sqrt(uu);
        for(uint32_t row=0;row<SV_cols;row++) {
            Q_t[col][row]=Q_t[col][row]/uu;
        }
    }
    for(uint32_t row=0;row<SV_rows;row++){
        for(uint32_t col=0;col<SV_cols;col++) {
            q[row][col]=Q_t[col][row];
        }
    }


    return q;
}

void qr(std::vector<std::vector<val_t>> const & mat,std::vector<std::vector<val_t>> & Q/*, std::vector<std::vector<val_t>> & R*/){
    int m = mat.size();
    int n = mat[0].size();

    // array of factor Q1, Q2, ... Qm
    std::vector<std::vector<std::vector<val_t>>> qv(m);

    // temp array
    std::vector<std::vector<val_t>> z;
    z=mat;
    std::vector<std::vector<val_t>> z1(m,std::vector<val_t>(n));

    for (int k = 0; k < n && k < m - 1; k++) {

        std::vector<val_t> e(m), x(m);
        val_t a;

        // compute minor
        z1=std::vector<std::vector<val_t>> (m,std::vector<val_t>(n));
        for (int i = 0; i < k; i++)
            z1[i][i] = 1.0;
        for (int i = k; i < z.size(); i++)
            for (int j = k; j < z[0].size(); j++)
                z1[i][j] = z[i][j];

        // extract k-th column into x
        a=0;
        for(uint32_t row=0;row<m;row++){
            val_t tmp;
            tmp = z1[row][k];
            x[row]= tmp;
            a+=tmp*tmp;
        }

        a = sqrt(a);
        if (mat[k][k] > 0) a = -a;

        for (int i = 0; i < e.size(); i++)
            e[i] = (i == k) ? 1 : 0;

        // e = x + a*e
        val_t span_e;
        span_e=0.0;
        for (int i = 0; i < e.size(); i++){
            val_t tmp;
            tmp = x[i]+a*e[i];
            e[i] = tmp;
            span_e+=tmp*tmp;
        }
        span_e=sqrt(span_e);

        // e = e / ||e||
        for (int i = 0; i < e.size(); i++)
            e[i]=e[i]/span_e;

        // qv[k] = I - 2 *e*e^T
        qv[k]=std::vector(e.size(),std::vector<val_t>(e.size()));
        for (int i = 0; i < e.size(); i++) {
            for (int j = 0; j < n; j++)
                qv[k][i][j] = -2 * e[i] * e[j];
        }
        for (int i = 0; i < e.size(); i++)
            qv[k][i][i] += 1;


        // z = qv[k] * z1
        z=std::vector(m,std::vector<val_t>(n));
        for (int ii = 0; ii < qv[k].size(); ii++)
            for (int jj = 0; jj < z1[0].size(); jj++)
                for (int kk = 0; kk < qv[k][0].size(); kk++)
                    z[ii][jj] += qv[k][ii][kk] * z1[kk][jj];

    }

    Q = qv[0];

    // after this loop, we will obtain Q (up to a transpose operation)
    for (int i = 1; i < n && i < m - 1; i++) {

        //z1.mult(qv[i], Q);
        z1=std::vector(m,std::vector<val_t>(n));
        for (int ii = 0; ii < qv[i].size(); ii++)
            for (int jj = 0; jj < Q[0].size(); jj++)
                for (int kk = 0; kk < qv[i][0].size(); kk++)
                    z1[ii][jj] += qv[i][ii][kk] * Q[kk][jj];

        Q = z1;

    }
    /*
    if(R){
        //R.mult(Q, mat);
        R=nullptr;
    }*/
    //Q.transpose();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < i; j++) {
            val_t t = Q[i][j];
            Q[i][j] = Q[j][i];
            Q[j][i] = t;
        }
    }

}
/*, std::vector<std::vector<val_t>> & R*/
/*void qrPP(std::vector<std::vector<val_t>> const & A,std::vector<std::vector<val_t>> & Q) {
    size_t dim;
    const size_t nrow = A.size();
    const size_t ncol = A[0].size();
    (ncol < nrow) ? dim = ncol : dim = nrow;

    val_t * Qi;
    Qi = new val_t[nrow*nrow];
    for(size_t i=0;i<nrow;i++){
        Qi[(i*nrow)+i]=0.0;
    }
    val_t * R;
    R = new val_t[nrow*ncol];
    for(size_t i=0;i<nrow;i++){
        for(size_t j=0;j<ncol;j++){
            R[i*ncol+j]=A[i][j];
        }
    }
    Matrix II;
    for (size_t i = 0; i != dim; i++)
    {
        val_t *  col;
        col = new val_t[A[0].size()];
        for(size_t k =0; k<A[0].size();k++)
            col[k] = A[0][k];
        val_t *  hous;
        hous = householder(col);
        val_t colNorm;
        colNorm=0.0;
        for(size_t k =0; k<A[0].size();k++)
            colNorm+=col[k]*col[k];
        colNorm=sqrt(colNorm);
        if (colNorm == 0.0)
        {
            col[0] = 1.0;
        }
        val_t *  e;
        e = new val_t[A[0].size()];
        for(size_t k =1; k<A[0].size();k++)
            e[k] = 0.0;
        e[0] = (col[0]/abs(col[0])) * colNorm;
        val_t *  v;
        v = new val_t[A[0].size()];
        for(size_t k =0; k<A[0].size();k++)
            v[k] = col[k] + e[k];
        val_t *  v;
        auto vT = transpose(v);
        auto vvT = v * vT;
        auto vTv = vT * v;
        auto H = identity_mat(dim) - ((2.0 * vvT) / vTv[0]);
        return H;



        if (i == 0)
            II = hous;
        else
        {
            II = identity_mat(nrow);
            II(range(i, nrow - 1), range(i, nrow - 1)) = hous;
        }
        *R = II * (*R);
        *Q = (*Q) * II;
        auto HA = hous * A;
        A = HA.sub_matrix(range(1, HA.get_nb_rows() - 1), range(1, HA.get_nb_columns() - 1));
    }
}*/
void qrP(std::vector<std::vector<val_t>> const & A,std::vector<std::vector<val_t>> & Q/*, std::vector<std::vector<val_t>> & R*/){
    // Make a copy of the input matrix.
    std::vector<std::vector<val_t>> inputMatrix;
    inputMatrix = A;

    // Verify that the input matrix is square.
    //if (A.size() != A[0].size())
        //return;//do nothing

    // Determine the number of columns (and rows, since the matrix is square).
    int numCols = A[0].size();
    int numRows = A.size();

    // Create a vector to store the P matrices for each column.
    std::vector<std::vector<std::vector<val_t>>> Plist;

    // Loop through each column.
    for (int j=0; j < numCols && j < numRows - 1; ++j)
    {
        // Create the a1 and b1 vectors.
        // a1 is the column vector from A.
        // b1 is the vector onto which we wish to reflect a1.
        std::vector<val_t> a1 (numRows-j);
        std::vector<val_t> b1 (numRows-j);
        for (int i=j; i<numRows; ++i)
        {
            a1[i-j]= inputMatrix[i][j];
            //b1[i-j]=0.0;
        }
        b1[0]=1.0;

        // Compute the norm of the a1 vector.
        val_t a1norm;
        a1norm=0;
        for(int ii=0;ii<a1.size();ii++)
            a1norm+= a1[ii]*a1[ii];
        a1norm=sqrt(a1norm);

        // Compute the sign we will use.
        int sgn = -1;
        if (a1[0] < 0.0)
            sgn = 1;

        // Compute the u-vector.
        std::vector<val_t> u (a1.size());
        val_t uNorm;
        uNorm=0;
        for(int ii=0;ii<a1.size();ii++){
            val_t tmp;
            tmp=a1[ii] - (sgn * a1norm * b1[ii]);
            u[ii] = tmp;
            uNorm+=tmp*tmp;
        }
        uNorm=sqrt(uNorm);

        // Compute the n-vector.
        std::vector<val_t> n (u.size());
        for(int ii=0;ii<u.size();ii++)
            n[ii] = u[ii]/uNorm;
        a1.clear();
        b1.clear();
        u.clear();
        // Compute Ptemp.
        std::vector<std::vector<val_t>> Ptemp;
        Ptemp=std::vector(numRows-j,std::vector<val_t>(numCols-j));
        for (int ii = 0; ii < numRows-j; ii++) {
            for (int jj = 0; jj < numCols-j; jj++)
                Ptemp[ii][jj] = -2 * n[ii] * n[jj];
        }
        n.clear();
        for (int ii = 0; ii < numRows-j && ii < numCols-j; ii++)
            Ptemp[ii][ii] += 1.0;

        // Form the P matrix with the original dimensions.
        std::vector<std::vector<val_t>> P (numRows, std::vector<val_t>(numCols));
        //P.SetToIdentity();
        for (int row=0; row<numRows && row < numCols; ++row)
            P[row][row]=1.0;
        for (int row=j; row<numRows; ++row)
        {
            for (int col=j; col<numCols; ++col)
            {
                P[row][col]=Ptemp[row-j][col-j];
            }
        }
        Ptemp.clear();

        // Store the result into the Plist vector.
        Plist.push_back(P);
        P.clear();

        // Apply this transform matrix to inputMatrix and use this result
        // next time through the loop.
        std::vector<std::vector<val_t>> rMat(Plist[j].size(),std::vector<val_t>(inputMatrix[0].size()));
        for (int lhsRow=0; lhsRow<Plist[j].size(); lhsRow++)
        {
            // Loop through each column on the RHS.
            for (int rhsCol=0; rhsCol<inputMatrix[0].size(); rhsCol++)
            {
                val_t elementResult;
                elementResult=0.0;
                // Loop through each element of this LHS row.
                for (int lhsCol=0; lhsCol<Plist[j][0].size(); lhsCol++)
                {
                    // Perform the calculation on these elements.
                    elementResult += Plist[j][lhsRow][lhsCol] * inputMatrix[lhsCol][rhsCol];
                }

                // Store the result.
                rMat[lhsRow][rhsCol] = elementResult;
            }
        }
        inputMatrix=rMat;
        rMat.clear();
    }

    // Compute Q.
    std::vector<std::vector<val_t>> Qmat;
    Qmat= Plist[0];
    for (int i=1; i<(numRows-1); ++i)
    {
        std::vector<std::vector<val_t>> tmpMat(Plist[i][0].size(),std::vector<val_t>(Plist[i].size()));
        for(int iRow=0; iRow<Plist[i][0].size(); iRow++){
            for(int iCol=0; iCol<Plist[i].size(); iCol++){
                tmpMat[iRow][iCol]=Plist[i][iCol][iRow];
            }
        }

        std::vector<std::vector<val_t>> rMat(Qmat.size(),std::vector<val_t>(tmpMat[0].size()));
        for (int lhsRow=0; lhsRow<Qmat.size(); lhsRow++)
        {
            // Loop through each column on the RHS.
            for (int rhsCol=0; rhsCol<tmpMat[0].size(); rhsCol++)
            {
                val_t elementResult;
                elementResult=0.0;
                // Loop through each element of this LHS row.
                for (int lhsCol=0; lhsCol<Qmat[0].size(); lhsCol++)
                {
                    // Perform the calculation on these elements.
                    elementResult += Qmat[lhsRow][lhsCol] * tmpMat[lhsCol][rhsCol];
                }

                // Store the result.
                rMat[lhsRow][rhsCol] = elementResult;
            }
        }
        Qmat=rMat;
        tmpMat.clear();
        rMat.clear();
    }

    // Return the Q matrix.
    Q = Qmat;
    Qmat.clear();
/*
    // Compute R.
    int numElements = Plist.size();
    qbMatrix2<T> Rmat = Plist.at(numElements-1);
    for (int i=(numElements-2); i>=0; --i)
    {
        Rmat = Rmat * Plist.at(i);
    }
    Rmat = Rmat * A;

    // And return the R matrix.
    R = Rmat;
*/
}

void HyperEF(spmv::io::CSCMatrix <val_t> ar, uint32_t L, uint32_t R){
    spmv::io::CSCMatrix <val_t> ar_new;
    spmv::io::CSCMatrix <val_t> idx_mat;


    std::vector<val_t> Neff(ar.num_rows);
    std::vector<val_t> W(ar.num_cols);
    std::fill(Neff.begin(), Neff.end(), 0);
    std::fill(W.begin(), W.end(), 1.0);

    for(uint32_t loop=L;loop<=L;loop++){
        uint32_t mx;
        uint32_t initial;
        uint32_t SmS;
        uint32_t interval;
        uint32_t Nrv;
        uint32_t Nsm;
        uint32_t Ntot;
        std::vector <val_t> Qvec;
        std::vector<std::vector<val_t>> Eratio;

        mx = ar.num_rows;

        spmv::io::CSCMatrix <val_t> A;
        //star expansion
        A = StarW(ar,W);
        //computing the smoothed vectors
        initial = 0;

        SmS = 300;

        interval = 20;

        Nrv = 1;

        Nsm = int(val_t(SmS - initial) / val_t(interval));

        Ntot = Nrv * Nsm;

        //Qvec = zeros(val_t64, 0);

        //Eratio = zeros(val_t64, length(ar), Ntot)

        //SV = zeros(val_t64, mx, Ntot)

        std::vector<std::vector<val_t>> SV(mx,std::vector<val_t>(Ntot));
        for(uint32_t ii=0;ii<Nrv;ii++){
            std::vector<std::vector<val_t>> sm;            //vector<vector<int>> vec( n , vector<int> (m, 0));
            //sm = zeros(mx, Nsm)
            std::vector<val_t> rv(A.num_rows);

            /*
            std::generate(rv.begin(), rv.end(), [&](){return ((val_t(get_random(0,INT32_MAX) / val_t(INT32_MAX)))- 0.5)*2.0;});
            */



            std::string line;
            std::ifstream myFile;            //creates stream myFile
            std::ofstream oFile;
            myFile.precision(19);
            oFile.precision(19);
            myFile.open("RVdata");  //opens .txt file
            oFile.open("RVdataO");
            uint32_t lines_read;
            lines_read=0;
            while (myFile.good()) {
                line.clear();
                getline(myFile, line);
                rv[lines_read++]=std::strtod(line.c_str(),nullptr);
                line=std::to_string(rv[lines_read-1]);
                oFile << rv[lines_read-1] << std::endl;
                if(lines_read==A.num_rows)
                    break;
            }
            myFile.close();
            oFile.close();


            sm = Filter(rv, SmS, A, mx, initial, interval, Nsm);

            std::ofstream outfile ("Debug");
            outfile.precision(15);
            for(uint32_t j=0;j<mx;j++){
                outfile << sm[j][Nsm-1] << std::endl;
            }
            outfile.close();

            //SV[:, (ii-1)*Nsm+1 : ii*Nsm] = sm
            for(uint32_t j=0;j<mx;j++){
                for(uint32_t i=0;i<Nsm;i++) {
                    SV[j][ii*Nsm+i] = sm[j][i];
                }
            }

        }
        std::vector<std::vector<val_t>> testA(3,std::vector<val_t>(3));
        testA[0]={12,-51,4};
        testA[1]={6,167,-68};
        testA[2]={-4,24,-41};
        testA[0]={0.5,0.75,0.5};
        testA[1]={1.0,0.5,0.75};
        testA[2]={0.25,0.25,0.25};
        std::vector<std::vector<val_t>> testQ;
        testQ=QR_decomposition_return_Q_only(testA);
        //https://github.com/QuantitativeBytes/qbLinAlg
        //qrP(testA,testQ);
        //std::vector<std::vector<val_t>> Q;
        qrP(SV,Q);
        std::ofstream out("cpp_filter.csv");
        out.precision(15);
        for (auto& row : Q) {
            for (auto col : row)
                out << col <<',';
            out << std::endl;
        }
        out.close();
    }
}
