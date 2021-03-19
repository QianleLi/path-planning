#include <hw_tool.h>

using namespace std;
using namespace Eigen;

void Homeworktool::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
}

void Homeworktool::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      
    
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool Homeworktool::isObsFree(const double coord_x, const double coord_y, const double coord_z)
{
    Vector3d pt;
    Vector3i idx;
    
    pt(0) = coord_x;
    pt(1) = coord_y;
    pt(2) = coord_z;
    idx = coord2gridIndex(pt);

    int idx_x = idx(0);
    int idx_y = idx(1);
    int idx_z = idx(2);

    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d Homeworktool::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i Homeworktool::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d Homeworktool::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

double Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position)
{
    double optimal_cost = 100000; // this just to initial the optimal_cost, you can delete it 
    /*
    STEP 2: go to the hw_tool.cpp and finish the function Homeworktool::OptimalBVP
    the solving process has been given in the document

    because the final point of trajectory is the start point of OBVP, so we input the pos,vel to the OBVP

    after finish Homeworktool::OptimalBVP, the Trajctory_Cost will record the optimal cost of this trajectory
    */
   //_start_position和_start_velocity是轨迹末端的位置和速度，那么轨迹开始的位置就是target_pt，而速度就是(0,0,0)
    double dPx = _target_position(0) - _start_position(0);  
    double dPy = _target_position(1) - _start_position(1);  
    double dPz = _target_position(2) - _start_position(2);  
    double dVx = _start_velocity(0);
    double dVy = _start_velocity(1);
    double dVz = _start_velocity(2);
    //计算多项式的系数u1T^4  + u2T^3+ u3T^2 + u4T + u5 = 0 = J'
    double u1 = 1; double u2 = 0;
    double u3 = 16 * (dVx*dVx + dVy*dVy + dVz*dVz);
    double u4 = -48 * (dPx*dVx + dPy*dVy + dPz*dVz);
    double u5 = 36 * (dPx*dPx + dPy*dPy + dPz*dPz);
    Eigen::Matrix<double, 4, 4>  companion = Eigen::Matrix<double, 4, 4>::Zero();
    companion (0,0) = -u2/u1;
    companion (0,1) = -u3/u1;
    companion (0,2) = -u4/u1;
    companion (0,3) = -u5/u1;
    companion (1,0) = 1.0;
    companion (2,1) = 1.0;
    companion (3,2) = 1.0;
    Eigen::EigenSolver<Eigen::Matrix<double, 4, 4>> eigen_solver (companion);
    //获取矩阵特征值 4*1
    MatrixXcd evals = eigen_solver.eigenvalues();
    //注意这里定义的MatrixXd里没有c
    MatrixXd evalsReal;
    //获取特征值实数部分
    evalsReal=evals.real();
    //特征值为两正两负，四个数绝对值一样
    //std::cout << evalsReal << std::endl;
    //std::cout << "========" << std::endl;
    double T;
    if(evalsReal(0) > evalsReal(2)) T = evalsReal(0);
    else T = evalsReal(2);
    double alpha1 = -12*dPx/(T*T*T) + 6*dVx/(T*T);
    double alpha2 = -12*dPy/(T*T*T) + 6*dVy/(T*T);
    double alpha3 = -12*dPy/(T*T*T) + 6*dVy/(T*T);
    double beta1 = 6*dPx/(T*T) - 2*dVx/T;
    double beta2 = 6*dPy/(T*T) - 2*dVy/T;
    double beta3 = 6*dPz/(T*T) - 2*dVz/T;
    optimal_cost = T + (alpha1*alpha1*T*T*T/3 + alpha1*beta1*T*T + beta1*beta1*T) + 
                                            (alpha2*alpha2*T*T*T/3 + alpha2*beta2*T*T + beta2*beta2*T) + 
                                            (alpha3*alpha3*T*T*T/3 + alpha3*beta3*T*T + beta3*beta3*T);
    return optimal_cost;
}