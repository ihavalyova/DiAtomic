#include <iostream>
#include <iomanip>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<double> calculate_wavenumbers_list(py::array_t<double> ulevels_inp, 
                                               py::array_t<double> llevels_inp)
{
    /* read input arrays buffer info */
    py::buffer_info ulevels_buf = ulevels_inp.request();
    py::buffer_info llevels_buf = llevels_inp.request();

    double *ulevesls_ptr = (double *) ulevels_buf.ptr;
    double *llevels_ptr = (double *) llevels_buf.ptr;

    int ushapex = ulevels_buf.shape[0];
    int ushapey = ulevels_buf.shape[1];
    int lshapex = llevels_buf.shape[0];
    int lshapey = llevels_buf.shape[1];
    int cw_shapex = ushapex * lshapex;
    int cw_shapey = 2*ushapey + 2;

    /* allocate the output buffer */
    int result_size = cw_shapex;

    int ncols = 18;  // number of columns in the output array
    py::array_t<double> result = py::array_t<double>(result_size*ncols);
    py::buffer_info result_buf = result.request();

    double *result_ptr = (double *) result_buf.ptr;

    /* allocate 2D arrays and fill them */
    double **ulevels = new double*[ushapex];
    for(int i = 0; i < ushapex; i++)
        ulevels[i] = new double[ushapey];

    for (int idx = 0; idx < ushapex; idx++)
        for (int idy = 0; idy < ushapey; idy++)
            ulevels[idx][idy] = ulevesls_ptr[idx*ushapey+idy];

    double **llevels = new double*[lshapex];
    for(int i = 0; i < lshapex; i++)
        llevels[i] = new double[lshapey];

    for (int idx = 0; idx < lshapex; idx++)
        for (int idy = 0; idy < lshapey; idy++)
            llevels[idx][idy] = llevels_ptr[idx*lshapey+idy];

    /* obtain the list of calculated wavenumbers */
    int wcount = 0;
    int nrows = 0;
    for(int i = 0; i < ushapex; i++)
    {
        auto ulevel = ulevels[i];

        int uind = (int)ulevel[0];
        int uv = (int)ulevel[1];
        double uj = ulevel[2];
        double uenrg = ulevel[3];      
        int up = (int)ulevel[4];
        int us = (int)ulevel[5];
        int ulam = ulevel[7];
        int uomg = (int)ulevel[8];

        for(int j = 0; j < lshapex; j++)
        {
            auto llevel = llevels[j];

            int lind = (int)llevel[0];
            int lv = (int)llevel[1];
            double lj = llevel[2];
            double lenrg = llevel[3];
            int lp = (int)llevel[4];
            int ls = (int)llevel[5];
            int llam = llevel[7];
            int lomg = (int)llevel[8];

            // -1 means the transition is not allowed by the strict 
            // selection rules for J and symmetry
            int branch = -1;

            // selection rule for Pee lines
            if (uj == lj-1 && up == 1 && lp == 1){
                branch = 0;
            }
            // selection rule for Pff lines
            if (uj == lj-1 && up == 0 && lp == 0){
                branch = 1;
            }
            // selection rule for Qef lines
            if (uj == lj && up == 1 && lp == 0){
                branch = 2;
            }
            // selection rule for Qfe lines
            if (uj == lj && up == 0 && lp == 1){
                branch = 3;
            }
            // selection rule for Ree lines
            if (uj == lj+1 && up == 1 && lp == 1){
                branch = 4;
            }
            // selection rule for Rff lines
            if (uj == lj+1 && up == 0 && lp == 0){
                branch = 5;
            }

            double cwaven = uenrg - lenrg;

            if(cwaven > 0.0)
            {
                result_ptr[wcount++] = uind;
                result_ptr[wcount++] = uv;
                result_ptr[wcount++] = uj;
                result_ptr[wcount++] = up;
                result_ptr[wcount++] = us;
                result_ptr[wcount++] = uenrg;
                result_ptr[wcount++] = ulam;
                result_ptr[wcount++] = uomg;
                result_ptr[wcount++] = lind;
                result_ptr[wcount++] = lv;
                result_ptr[wcount++] = lj;
                result_ptr[wcount++] = lp;
                result_ptr[wcount++] = ls;
                result_ptr[wcount++] = lenrg;
                result_ptr[wcount++] = llam;
                result_ptr[wcount++] = lomg;
                result_ptr[wcount++] = cwaven;
                result_ptr[wcount++] = branch;
                nrows++;
            }
        }
    }

    /* reshape result */
    result.resize({nrows, ncols});

    /* delete allocated memory */
    for(int i = 0; i < ushapex; i++)
        delete [] ulevels[i];
    delete [] ulevels;

    for(int i = 0; i < lshapex; i++)
        delete [] llevels[i];
    delete [] llevels;

    return result;
}

py::array_t<double> calculate_wavenumbers_list_with_obs(py::array_t<double, py::array::c_style> ulevels_inp, 
                                                        py::array_t<double, py::array::c_style> llevels_inp,
                                                        py::array_t<double, py::array::c_style> ewaven_inp)
{
    /* read input arrays buffer info */
    py::buffer_info ulevels_buf = ulevels_inp.request();
    py::buffer_info llevels_buf = llevels_inp.request();
    py::buffer_info ewaven_buf = ewaven_inp.request();

    double *ulevels_ptr = (double *) ulevels_buf.ptr;
    double *llevels_ptr = (double *) llevels_buf.ptr;
    double *ewaven_ptr = (double *) ewaven_buf.ptr;

    int ushapex = ulevels_buf.shape[0];
    int ushapey = ulevels_buf.shape[1];
    int lshapex = llevels_buf.shape[0];
    int lshapey = llevels_buf.shape[1];
    int ew_shapex = ewaven_buf.shape[0];
    int ew_shapey = ewaven_buf.shape[1];
    int cw_shapex = ushapex * lshapex;
    int cw_shapey = 17; // 2*ushapey + 2;

    /* allocate the output buffer */
    int result_size = 0;
    if (ew_shapex > cw_shapex)
        result_size = ew_shapex;
    else
        result_size = cw_shapex;

    int ncols = 23;  // number of columns in the output array
    py::array_t<double> result = py::array_t<double>(result_size*ncols);
    py::buffer_info result_buf = result.request();

    double *result_ptr = (double *) result_buf.ptr;

    /* allocate 2D arrays and fill them */
    double **ulevels = new double*[ushapex];
    for(int i = 0; i < ushapex; i++)
        ulevels[i] = new double[ushapey];

    for (int idx = 0; idx < ushapex; idx++)
        for (int idy = 0; idy < ushapey; idy++)
            ulevels[idx][idy] = ulevels_ptr[idx*ushapey+idy];

    double **llevels = new double*[lshapex];
    for(int i = 0; i < lshapex; i++)
        llevels[i] = new double[lshapey];

    for (int idx = 0; idx < lshapex; idx++)
        for (int idy = 0; idy < lshapey; idy++)
            llevels[idx][idy] = llevels_ptr[idx*lshapey+idy];

    double **waven_cal = new double*[cw_shapex];
    for(int i = 0; i < cw_shapex; i++)
        waven_cal[i] = new double[cw_shapey];

    double **waven_exp = new double*[ew_shapex];
    for(int i = 0; i < ew_shapex; i++)
        waven_exp[i] = new double[ew_shapey];

    for (int idx = 0; idx < ew_shapex; idx++)
        for (int idy = 0; idy < ew_shapey; idy++)
            waven_exp[idx][idy] = ewaven_ptr[idx*ew_shapey+idy];

    /* obtain the list of calculated wavenumbers */
    int wcount = 0;
    for(int i = 0; i < ushapex; i++)
    {
        auto ulevel = ulevels[i];

        int uind = (int)ulevel[0];
        int uv = (int)ulevel[1];
        double uj = ulevel[2];
        double uenrg = ulevel[3];      
        int up = (int)ulevel[4];
        int us = (int)ulevel[5];
        int ulam = ulevel[7];
        int uomg = (int)ulevel[8];

        for(int j = 0; j < lshapex; j++)
        {
            auto llevel = llevels[j];

            int lind = (int)llevel[0];
            int lv = (int)llevel[1];
            double lj = llevel[2];
            double lenrg = llevel[3];
            int lp = (int)llevel[4];
            int ls = (int)llevel[5];
            int llam = llevel[7];
            int lomg = (int)llevel[8];

            double cwaven = uenrg - lenrg;

            waven_cal[wcount][0] = uind;
            waven_cal[wcount][1] = uv;
            waven_cal[wcount][2] = uj;
            waven_cal[wcount][3] = up;
            waven_cal[wcount][4] = us;
            waven_cal[wcount][5] = uenrg;
            waven_cal[wcount][6] = ulam;
            waven_cal[wcount][7] = uomg;
            waven_cal[wcount][8] = lind;
            waven_cal[wcount][9] = lv;
            waven_cal[wcount][10] = lj;
            waven_cal[wcount][11] = lp;
            waven_cal[wcount][12] = ls;
            waven_cal[wcount][13] = lenrg;
            waven_cal[wcount][14] = llam;
            waven_cal[wcount][15] = lomg;
            waven_cal[wcount][16] = cwaven;
            wcount++;
        }
    }

    /* obtain the final list of wavenumbers */
    int count = 0;
    int nrows = 0;
    for(int i = 0; i < cw_shapex; i++)
    {
        auto cline = waven_cal[i];

        int cuind = (int)cline[0];
        int cuv = (int)cline[1];
        double cuj = cline[2];
        int cup = (int)cline[3];
        int cus = (int)cline[4];
        double cue = cline[5];
        int culam = cline[6];
        double cuomg = cline[7];

        int clind = (int)cline[8];
        int clv = (int)cline[9];
        double clj = cline[10];
        int clp = (int)cline[11];
        int cls = (int)cline[12];
        double cle = cline[13];
        int cllam = cline[14];
        double clomg = cline[15];
        double cwaven = cline[16];

        for(int j = 0; j < ew_shapex; j++)
        {
            auto eline = waven_exp[j];

            int euv = (int)eline[0];
            double euj = eline[1];
            int eup = (int)eline[2]; 
            int eus = (int)eline[3];
            int elv = (int)eline[4];
            double elj = eline[5];
            int elp = (int)eline[6];
            int els = (int)eline[7];
            double ewaven = eline[8];
            double eunc_wavens = eline[9];
            double eintens = eline[10];
            double eunc_intens = eline[11];

            if (cuv == euv && clv == elv && cuj == euj && clj == elj && 
                cup == eup && clp == elp && cus == eus && cls == els)
            {
                // -1 means the transition is not allowed by the strict 
                // selection rules for J and symmetry
                int branch = -1;

                // selection rule for Pee lines
                if (euj == elj-1 && eup == 1 && elp == 1){
                    branch = 0;
                }
                // selection rule for Pff lines
                if (euj == elj-1 && eup == 0 && elp == 0){
                    branch = 1;
                }
                // selection rule for Qef lines
                if (euj == elj && eup == 1 && elp == 0){
                    branch = 2;
                }
                // selection rule for Qfe lines
                if (euj == elj && eup == 0 && elp == 1){
                    branch = 3;
                }
                // selection rule for Ree lines
                if (euj == elj+1 && eup == 1 && elp == 1){
                    branch = 4;
                }
                // selection rule for Rff lines
                if (euj == elj+1 && eup == 0 && elp == 0){
                    branch = 5;
                }

                double diff_waven = ewaven - cwaven;

                result_ptr[count++] = cuind;
                result_ptr[count++] = euv;
                result_ptr[count++] = euj;
                result_ptr[count++] = eup;
                result_ptr[count++] = eus;
                result_ptr[count++] = cue;
                result_ptr[count++] = culam;
                result_ptr[count++] = cuomg;
                result_ptr[count++] = clind;
                result_ptr[count++] = elv;
                result_ptr[count++] = elj;
                result_ptr[count++] = elp;
                result_ptr[count++] = els;
                result_ptr[count++] = cle;
                result_ptr[count++] = cllam;
                result_ptr[count++] = clomg;
                result_ptr[count++] = cwaven;
                result_ptr[count++] = ewaven;
                result_ptr[count++] = diff_waven;
                result_ptr[count++] = eunc_wavens;
                result_ptr[count++] = eintens;
                result_ptr[count++] = eunc_intens;
                result_ptr[count++] = branch;
                nrows++;
            }
        }
    }

    /* reshape result */
    result.resize({nrows, ncols});

    /* delete allocated memory */
    for(int i = 0; i < ushapex; i++)
        delete [] ulevels[i];
    delete [] ulevels;

    for(int i = 0; i < lshapex; i++)
        delete [] llevels[i];
    delete [] llevels;

    for(int i = 0; i < cw_shapex; i++)
        delete [] waven_cal[i];
    delete [] waven_cal;

    for(int i = 0; i < ew_shapex; i++)
        delete [] waven_exp[i];
    delete [] waven_exp;

    return result;
}

PYBIND11_MODULE(wavenumbers, m)
{
    m.def("calculate_wavenumbers_list", 
          &calculate_wavenumbers_list, 
          "Calculate the final wavenumbers list");
    
    m.def("calculate_wavenumbers_list_with_obs", 
          &calculate_wavenumbers_list_with_obs, 
          "Calculate the final wavenumbers list with observations");
}
