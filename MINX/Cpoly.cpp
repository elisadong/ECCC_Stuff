#include <Python.h>
#include <cmath>
#include <limits>
#include <iostream>
// #include <algorithm>

// did this NEED to be done in full cpp instead of like cython or something? No. Was it anyways? Yes. I moved on, so should you.

// TODO needs to work with just two passed values, maybe get center in this guy or something idk
// until its fixed, sorting is offically NotImplementedError. (my favourite kind of error)
// bool sort_by_atan2 (double px1, double py1, double cx1, double cy1, double px2, double py2, double cx2, double cy2){
//     return std::atan2(py1 - cy1, px1 - cx1) > std::atan2(py2 - cy2, px2 - cx2);
// }


static PyObject *polygon_check_single(PyObject *self, PyObject *args, PyObject *kwargs){

    // std::cout << "--STARTING SINGLE CHECK\n\n";
    // std::cout << "----PARSING ARGS\n";

    Py_Initialize();

    PyObject *listObj, *tupObj, *pntObj;

    int *ccwsorted = NULL;

    static char *kwlist[] = {"polypoints", "point", "ccwsorted", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|p", kwlist, &PyList_Type, &listObj, &PyTuple_Type, &pntObj, &ccwsorted)) {
        return NULL;
    }

    // std::cout << "------TUPLE AND KWARGS READ\n";

    double point[2];
    point[0] = (double)  PyFloat_AsDouble(PyTuple_GetItem(pntObj, 0));
    point[1] = (double) PyFloat_AsDouble(PyTuple_GetItem(pntObj, 1));

    int N = PyList_Size(listObj);

    // std::cout << "------STARTING POLYPOINTS INIT\n";
    // std::cout << "--------N: " << N << "\n";
    // std::cout << "--------ccwsorted: " << ccwsorted << "\n";
    // std::cout << "--------point: (" << point[0] << ' ' << point[1] << ")\n";
    double** polypoints = new double*[N];
    for(int i=0; i < N; ++i){
        polypoints[i] = new double[2];

        tupObj = PyList_GetItem(listObj, i);
        polypoints[i][0] = (double) PyFloat_AsDouble(PyTuple_GetItem(tupObj, 0));
        polypoints[i][1] = (double) PyFloat_AsDouble(PyTuple_GetItem(tupObj, 1));
        // std::cout << "--------poly: (" << polypoints[i][0] << ' ' << polypoints[i][1] << ")\n";
    };
    // std::cout << "------ENDING POLYPOINTS INIT\n";

   
    // std::cout << "----ENDING ARG PARSING\n\n";
    /**************************************************************************************/

    // std::cout << "----STARTING SORT AND OBV CHECK\n";

    // (1) sort and obv check

    // init top, bot, left, right to just the first elem, will check all anyways
    double left, right, top, bot;
    left = right = polypoints[0][0];
    top = bot = polypoints[0][1];

    for (int i=0; i < N; ++i){
        if (polypoints[i][0] == point[0] && polypoints[i][1] == point[1]) { Py_RETURN_TRUE; };

        if (polypoints[i][0] < left) { left = polypoints[i][0]; };
        if (polypoints[i][0] > right) { right = polypoints[i][0]; };
        if (polypoints[i][1] < bot) { bot = polypoints[i][1]; };
        if (polypoints[i][1] > top) { top = polypoints[i][1]; };
    };

    // std::cout << "------top: (" << top << ")\n";
    // std::cout << "------bot: (" << bot << ")\n";
    // std::cout << "------left: (" << left << ")\n";
    // std::cout << "------right: (" << right << ")\n";

    if (!(left < point[0] && point[0] < right) || !(bot < point[1] && point[1] < top)){
        // std::cout << "~~OBV CHECK FALSE -> RETURNING FALSE~~\n";
        Py_RETURN_FALSE;
    };
    // std::cout << "----ENDING SORT AND OBV CHECK\n\n";

    // std::cout << "----STARTING CCWSORT\n";

    // (2) ccwsort (DOESNT WORK YET DONT TRY IT PLZ)(HONESTLY YOU SHOULD ALWAYS DO YOUR OWN SORTING, DONT TRUST ME WITH THIS)

    if (!ccwsorted) {
        // std::cout << "------HASNT BEEN SORTED (ERROR)\n";
        double cent_x = 0.0;  // Cant get std sum to work right so whatever
        double cent_y = 0.0;
        for (int i=0; i < N; ++i){
            cent_x += polypoints[i][0];
            cent_y += polypoints[i][1];
        }
        cent_x /= N;
        cent_y /= N;

        // std::sort(polypoints, polypoints + N, sort_by_atan2);  // TODO

    }
    // std::cout << "----ENDING CCWSORT\n\n";

    // (3) draw lines
    // std::cout << "----STARTING DRAW LINE\n";

    double slopes[N], intercepts[N];
    double xranges[N][2], yranges[N][2];

    // assign values instead of calling, and init with last point to avoid negative indexing (which is apparently a no go in C)
    // draw lines from [-1]->[0], [0]->[1], etc

    double x1 = polypoints[N-1][0];
    double y1 = polypoints[N-1][1];

    for (int i=0; i < N; ++i){
        double x2 = polypoints[i][0];
        double y2 = polypoints[i][1];


        // std::cout << "------line " << i << " : ";
        if (x2 - x1){
            slopes[i] = (y2 - y1) / (x2 - x1);
            intercepts[i] = y2 - ((y2 - y1) / (x2 - x1)) * x2;

            yranges[i][0] = std::numeric_limits<double>::quiet_NaN();
        } else {
            slopes[i] = std::numeric_limits<double>::quiet_NaN();
            intercepts[i] = std::numeric_limits<double>::quiet_NaN();
            
            // janky sorting cause std::sort keeps breaking
            yranges[i][0] = (y2 < y1) ? y2 : y1;
            yranges[i][1] = (y2 < y1) ? y1 : y2;
        }

        xranges[i][0] = (x2 < x1) ? x2 : x1;
        xranges[i][1] = (x2 < x1) ? x1 : x2;

        x1 = x2; y1 = y2;

        // std::cout << slopes[i] << "x+" << intercepts[i] << " with ranges (x,y): (" \
            << xranges[i][0] << ", " << xranges[i][1] << "), (" << yranges[i][0] << ", " << yranges[i][1] << ")\n";
    }
    // std::cout << "----ENDING DRAW LINE\n\n";

    // (4) check intercept
    // std::cout << "----STARTING INT CHECK\n";

    int results[N];

    for (int i=0; i < N; ++i){
        // vertical line checking
        if (!std::isnan(yranges[i][0])){
            if (point[0] == xranges[i][0] && yranges[i][0] < point[1] && point[1] < yranges[i][1]){
                // std::cout << "~~VERT CHECK TRUE -> RETURNING TRUE~~\n";
                Py_RETURN_TRUE;
            } else {
                // out of bounds result
                results[i] = -1;
            }
        // calc line for final check
        } else if (xranges[i][0] <= point[0] && point[0] <= xranges[i][1]){
            double line = slopes[i] * point[0] + intercepts[i];
            results[i] = line > point[1];
        }
        else {
            // out of bounds result
            results[i] = -1;
        }
        // std::cout << "------result " << i << " : " << results[i] << "\n";
    }
    // std::cout << "----ENDING INT CHECK\n\n";

    // (5) final check

    // std::cout << "----STARTING FINAL CHECK\n";

    bool allFalse = true;
    bool allTrue = true;

    for (int i=0; i < N; ++i){
        // skip the out of bound results
        if (results[i] < 0) continue;
        
        if (results[i]){
            allFalse = false;
        } else if (!results[i]){
            allTrue = false;
        }
    }
    // std::cout << "------allFalse: " << allFalse << "\n";
    // std::cout << "------allTrue: " << allTrue << "\n";

    for (int i=0; i < N; ++i){ delete[] polypoints[i]; }
    delete[] polypoints;

    if (allTrue || allFalse){
        // std::cout << "~~FINAL CHECK FALSE -> RETURNING FALSE~~\n";
        Py_RETURN_FALSE;
    } else {
        // std::cout << "~~FINAL CHECK TRUE -> RETURNING TRUE~~\n";
        Py_RETURN_TRUE;
    }

    /**************************************************************************************/
};


static PyObject *polygon_check(PyObject *self, PyObject *args, PyObject *kwargs){

    // std::cout << "--STARTING CHECK\n\n";
    // std::cout << "----PARSING ARGS\n";

    Py_Initialize();

    PyObject *polylistObj, *boxlistObj, *tupObj;

    int *ccwsorted = NULL;

    static char *kwlist[] = {"polypoints", "point", "ccwsorted", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|p", kwlist, &PyList_Type, &polylistObj, &PyList_Type, &boxlistObj, &ccwsorted)) {
        return NULL;
    }

    int poly_N = PyList_Size(polylistObj);
    int box_N = PyList_Size(boxlistObj);

    // std::cout << "------TUPLE AND KWARGS READ\n";
    // std::cout << "--------ccwsorted: " << ccwsorted << "\n";
    // std::cout << "--------poly_N: " << poly_N << "\n";
    // std::cout << "--------box_N: " << box_N << "\n";


    // std::cout << "------STARTING POLYPOINTS INIT\n";
    double** polypoints = new double*[poly_N];
    for(int i=0; i < poly_N; ++i){
        polypoints[i] = new double[2];

        tupObj = PyList_GetItem(polylistObj, i);
        polypoints[i][0] = (double) PyFloat_AsDouble(PyTuple_GetItem(tupObj, 0));
        polypoints[i][1] = (double) PyFloat_AsDouble(PyTuple_GetItem(tupObj, 1));
        // std::cout << "--------poly: (" << polypoints[i][0] << ' ' << polypoints[i][1] << ")\n";
    };
    // std::cout << "------ENDING POLYPOINTS INIT\n";

    // std::cout << "------STARTING BBOX INIT\n";
    double** bbox = new double*[box_N];
    for(int i=0; i < box_N; ++i){
        bbox[i] = new double[2];

        tupObj = PyList_GetItem(boxlistObj, i);
        bbox[i][0] = (double) PyFloat_AsDouble(PyTuple_GetItem(tupObj, 0));
        bbox[i][1] = (double) PyFloat_AsDouble(PyTuple_GetItem(tupObj, 1));
        // std::cout << "--------box: (" << bbox[i][0] << ' ' << bbox[i][1] << ")\n";
    };
    // std::cout << "------ENDING BBOX INIT\n";

   
    // std::cout << "----ENDING ARG PARSING\n\n";
    /**************************************************************************************/

    // std::cout << "----STARTING SORT AND OBV CHECK\n";

    // (1) sort and obv check

    for (int i=1; i < poly_N; ++i){
        for (int j=1; j < box_N; ++j){
            if (polypoints[i][0] == bbox[j][0] && polypoints[i][1] == bbox[j][1]) { 
                // std::cout << "~~POLYGONS SHARE CORNER -> RETURNING TRUE~~\n";
                Py_RETURN_TRUE;
            };
        };
    };

    // init top, bot, left, right to just the first elem, will check all anyways
    double p_left, p_right, p_top, p_bot, b_left, b_right, b_top, b_bot;

    p_left = p_right = polypoints[0][0];
    p_top = p_bot = polypoints[0][1];

    b_left = b_right = bbox[0][0];
    b_top = b_bot = bbox[0][1];

    for (int i=1; i < poly_N; ++i){
        if (polypoints[i][0] < p_left) { p_left = polypoints[i][0]; };
        if (polypoints[i][0] > p_right) { p_right = polypoints[i][0]; };
        if (polypoints[i][1] < p_bot) { p_bot = polypoints[i][1]; };
        if (polypoints[i][1] > p_top) { p_top = polypoints[i][1]; };
    };
    for (int i=1; i < box_N; ++i){
        if (bbox[i][0] < b_left) { b_left = bbox[i][0]; };
        if (bbox[i][0] > b_right) { b_right = bbox[i][0]; };
        if (bbox[i][1] < b_bot) { b_bot = bbox[i][1]; };
        if (bbox[i][1] > b_top) { b_top = bbox[i][1]; };
    };

    // std::cout << "------p_top: (" << p_top << ")\n";
    // std::cout << "------p_bot: (" << p_bot << ")\n";
    // std::cout << "------p_left: (" << p_left << ")\n";
    // std::cout << "------p_right: (" << p_right << ")\n";

    // std::cout << "------b_top: (" << b_top << ")\n";
    // std::cout << "------b_bot: (" << b_bot << ")\n";
    // std::cout << "------b_left: (" << b_left << ")\n";
    // std::cout << "------b_right: (" << b_right << ")\n";

    if (b_top < p_bot || b_bot > p_top || b_right < p_left || b_left > p_right){
        // std::cout << "~~BBOX OUTSIDE POLY CROP -> RETURNING FALSE~~\n";
        Py_RETURN_FALSE;
    } else if (b_top > p_top && b_bot < p_bot && b_right > p_right && b_left < p_left){
        // std::cout << "~~BBOX SURROUNDS POLY -> RETURNING TRUE~~\n";
        Py_RETURN_TRUE;
    }

    // std::cout << "----ENDING SORT AND OBV CHECK\n\n";

    // std::cout << "----STARTING CCWSORT\n";

    // (2) ccwsort (DOESNT WORK YET DONT TRY IT PLZ)(HONESTLY YOU SHOULD ALWAYS DO YOUR OWN SORTING, DONT TRUST ME WITH THIS)

    if (!ccwsorted) {
        // std::cout << "------HASNT BEEN SORTED (ERROR)\n";
    }
    // std::cout << "----ENDING CCWSORT\n\n";

    // (3) draw lines
    // std::cout << "----STARTING DRAW LINE\n";

    double p_slopes[poly_N], p_intercepts[poly_N];
    double p_xranges[poly_N][2], p_yranges[poly_N][2];

    // assign values instead of calling, and init with last point to avoid negative indexing (which is apparently a no go in C)
    // draw lines from [-1]->[0], [0]->[1], etc

    double px1 = polypoints[poly_N-1][0];
    double py1 = polypoints[poly_N-1][1];

    for (int i=0; i < poly_N; ++i){
        double px2 = polypoints[i][0];
        double py2 = polypoints[i][1];

        // std::cout << "------p_line ("<< i-1 << '-' << i << ") : ";
        if (px2 - px1){
            p_slopes[i] = (py2 - py1) / (px2 - px1);
            p_intercepts[i] = py2 - ((py2 - py1) / (px2 - px1)) * px2;

            p_yranges[i][0] = std::numeric_limits<double>::quiet_NaN();
        } else {
            p_slopes[i] = std::numeric_limits<double>::quiet_NaN();
            p_intercepts[i] = std::numeric_limits<double>::quiet_NaN();
            
            // janky sorting cause std::sort keeps breaking
            p_yranges[i][0] = (py2 < py1) ? py2 : py1;
            p_yranges[i][1] = (py2 < py1) ? py1 : py2;
        }

        p_xranges[i][0] = (px2 < px1) ? px2 : px1;
        p_xranges[i][1] = (px2 < px1) ? px1 : px2;

        px1 = px2; py1 = py2;

        // std::cout << p_slopes[i] << "x+" << p_intercepts[i] << " with ranges (x,y): (" \
            << p_xranges[i][0] << ", " << p_xranges[i][1] << "), (" << p_yranges[i][0] << ", " << p_yranges[i][1] << ")\n";
    };

    double b_slopes[box_N], b_intercepts[box_N];
    double b_xranges[box_N][2], b_yranges[box_N][2];

    // assign values instead of calling, and init with last point to avoid negative indexing (which is apparently a no go in C)
    // draw lines from [-1]->[0], [0]->[1], etc

    double bx1 = bbox[box_N-1][0];
    double by1 = bbox[box_N-1][1];

    for (int i=0; i < box_N; ++i){
        double bx2 = bbox[i][0];
        double by2 = bbox[i][1];

        // std::cout << "------b_line ("<< i-1 << '-' << i << ") : ";
        if (bx2 - bx1){
            b_slopes[i] = (by2 - by1) / (bx2 - bx1);
            b_intercepts[i] = by2 - ((by2 - by1) / (bx2 - bx1)) * bx2;

            b_yranges[i][0] = std::numeric_limits<double>::quiet_NaN();
        } else {
            b_slopes[i] = std::numeric_limits<double>::quiet_NaN();
            b_intercepts[i] = std::numeric_limits<double>::quiet_NaN();
            
            // janky sorting cause std::sort keeps breaking
            b_yranges[i][0] = (by2 < by1) ? by2 : by1;
            b_yranges[i][1] = (by2 < by1) ? by1 : by2;
        }

        b_xranges[i][0] = (bx2 < bx1) ? bx2 : bx1;
        b_xranges[i][1] = (bx2 < bx1) ? bx1 : bx2;

        bx1 = bx2; by1 = by2;

        // std::cout << b_slopes[i] << "x+" << b_intercepts[i] << " with ranges (x,y): (" \
            << b_xranges[i][0] << ", " << b_xranges[i][1] << "), (" << b_yranges[i][0] << ", " << b_yranges[i][1] << ")\n";
    };

    // std::cout << "----ENDING DRAW LINE\n\n";

    // std::cout << "----STARTING INT CHECK\n";

    for (int i=0; i < poly_N; ++i){
        for (int j=0; j < box_N; ++j){

            // std::cout << "\n------line (p#, b#, pline, bline): ("<< i-1 << '-' << i << "), (" << j-1 << '-' << j << "), (" \
                        << p_slopes[i] << "x+" << p_intercepts[i] << "), (" << b_slopes[i] << "x+" << b_intercepts[i] << ")\n";

            if ((p_slopes[i] == b_slopes[j]) || (std::isnan(p_slopes[i]) && std::isnan(b_slopes[j]))){
                // std::cout << "--------lines are parallel\n";

                if (std::isnan(p_slopes[i]) && p_xranges[i][0] == b_xranges[j][0]){
                    // std::cout << "~~TWO LINES VERTICAL WITH SAME X -> RETURNING TRUE~~\n";
                    Py_RETURN_TRUE;  // both vertical, with same x
                } else { continue; }

            } else {
                // std::cout << "--------lines are *not* parallel\n";

                if (!std::isnan(p_yranges[i][0])){
                    // verical poly lines
                    double x_int = p_xranges[i][0];
                    double int_line = b_slopes[j] * x_int + b_intercepts[j];
                    // std::cout << "--------x int: " << x_int << "\n";
                    // std::cout << "--------int_line: " << int_line << "\n";
                    // std::cout << "--------p_yranges: (" << p_yranges[i][0] << ", " << p_yranges[i][1] << ")\n";
                    if ((p_yranges[i][0] <= int_line && int_line <= p_yranges[i][1]) && (b_xranges[j][0] <= x_int && x_int <= b_xranges[j][1])){
                        // std::cout << "~~VERTICAL POLY LINE IS INT -> RETURNING TRUE~~\n";
                        Py_RETURN_TRUE;
                    }

                } else if (!std::isnan(b_yranges[i][0])){
                    // vertical box lines
                    double x_int = b_xranges[i][0];
                    double int_line = p_slopes[i] * x_int + p_intercepts[i];
                    // std::cout << "--------x int: " << x_int << "\n";
                    // std::cout << "--------int_line: " << int_line << "\n";
                    // std::cout << "--------b_yranges: (" << b_yranges[j][0] << ", " << b_yranges[j][1] << ")\n";
                    if ((b_yranges[j][0] <= int_line && int_line <= b_yranges[j][1]) && (p_xranges[j][0] <= x_int && x_int <= p_xranges[j][1])){
                        // std::cout << "~~VERTICAL BOX LINE IS INT -> RETURNING TRUE~~\n";
                        Py_RETURN_TRUE;
                    }

                } else {
                    // std::cout << "--------neither line vertical, checking int\n";
                    // non-vertical lines (usual case)
                    double x_int = (b_intercepts[j] - p_intercepts[i]) / (p_slopes[i] - b_slopes[j]);

                    // std::cout << "--------x int: " << x_int << "\n";
                    // std::cout << "--------p_xrange: " << p_xranges[i][0] << '-' << p_xranges[i][1] << "\n";
                    // std::cout << "--------b_xrange: " << b_xranges[j][0] << '-' << b_xranges[j][1] << "\n";

                    if ((p_xranges[i][0] <= x_int && x_int <= p_xranges[i][1]) && (b_xranges[j][0] <= x_int && x_int <= b_xranges[j][1])){
                        // std::cout << "~~LINES INTERSECT -> RETURNING TRUE~~\n";
                        Py_RETURN_TRUE;
                    }
                }
            }
        };
    };
    // std::cout << "----ENDING INT CHECK\n";

    // std::cout << "--ENDING CHECK\n\n";

    // std::cout << "--CALLING SINGLE CHECK\n";

    PyObject *point = Py_BuildValue("(dd)", bbox[0][0], bbox[0][1]);
    // std::cout << "----point: " << point << "\n";
    PyObject *newargs = Py_BuildValue("OO", polylistObj, point);

    polygon_check_single(self, newargs, kwargs);
};


/*PyDoc_STRVAR(
    single_doc,
    "Function for checking if a single point lies within a given N-sided\n\
    convex polygon. (Concave polygons not possible without prior info,\n\
    please sort concave polygons before passing points to function)\n\n\
    WARNING: Cpoly does not currently support sorting, please pre-sort polygons.\n\
    Input: list of N tuples (lon, lat) of polygon corners\n\
           Tuple (lon, lat) of point to check\n\
           Bool, wether polygon has been pre-sorted counterclockwise, from pi rad\n\
    Output: Bool True if point lies within polygon else False"
);
*/


static PyMethodDef Cpoly_methods[] = {
    {"polygon_check_single", polygon_check_single, METH_VARARGS | METH_KEYWORDS, "single docstring..."},
    {"polygon_check", polygon_check, METH_VARARGS | METH_KEYWORDS, "bbox docstring..."},
    {}        /* Sentinel */
};

static struct PyModuleDef Cpoly = {
    PyModuleDef_HEAD_INIT,
    "Cpoly",   /* name of module */
    "mod doc",      /* module documentation, may be NULL */
    -1,        /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    Cpoly_methods
};

PyMODINIT_FUNC PyInit_Cpoly(void){
    Py_Initialize();
    return PyModule_Create(&Cpoly);
};
