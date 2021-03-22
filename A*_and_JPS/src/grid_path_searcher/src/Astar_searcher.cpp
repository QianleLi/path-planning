#include "../include/Astar_searcher.h"

using namespace std;
using namespace Eigen;

void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   //初始化
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    //id的最大值，用来进行边界检测
    GLX_SIZE = max_x_id;   
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;   
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    //用边界来初始化这个三维指针数组，代表着三维空间
    //因此所有的节点都在这里创建好了，在寻找邻点的函数AstarGetSucc里不能用新建的GridNodePtr
    //否则将无法判定哪个节点已经在closed或open list中
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids()
{   
    for(int i=0; i < GLX_SIZE ; i++)
        for(int j=0; j < GLY_SIZE ; j++)
            for(int k=0; k < GLZ_SIZE ; k++)
                resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

vector<Vector3d> AstarPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++){   
                //if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and close list
                if(GridNodeMap[i][j][k]->id == -1)  // visualize nodes in close list only
                    visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i AstarPathFinder::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i & index) const
{
    return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return  (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
            (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool AstarPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{  //边界检测，判断是不是障碍物
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets)
{   
    neighborPtrSets.clear();
    edgeCostSets.clear();
    /*
    *
    STEP 4: finish AstarPathFinder::AstarGetSucc yourself 
    please write your code below
    *
    *
    */
   //创建三个vector，用这三个vector来获得当前点的邻点，这里增加的是index的值
    const vector<int> deltaX{-1, -1, -1,  0,   0,  0,  -1,  1,  1,  1,  1, -1,  0,  1, -1,  0,  1,   0, -1,   0,  1,   1,   1,  0, -1, -1};
    const vector<int> deltaY{  0,   0, -1,  0, -1, -1, -1,  0, -1, -1, 0,   1,  1,  1,   1,  1,  1,  0, -1, -1,  -1,  0,   1,  1,   1,   0};
    const vector<int> deltaZ{  0,  1,   1,  1,   1,  0,   0,  0,  0,   1, 1,   1,  1,  1,   0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    for(size_t i=0; i<deltaX.size(); i++)
    {
        Vector3i index(currentPtr->index(0)+deltaX[i], currentPtr->index(1)+deltaY[i], currentPtr->index(2)+deltaZ[i]);
        
        if(isFree(index))
        {
            //这里要用GridNodeMap中已经存在的指针，而不能创建新的
            
            double new_cost = sqrt(pow(deltaX[i],2) + pow(deltaY[i],2) + pow(deltaZ[i],2));
            neighborPtrSets.push_back(GridNodeMap[index(0)][index(1)][index(2)]);
            //new_cost是用coord计算出来的，要用resolution来转换到index
            edgeCostSets.push_back(new_cost * resolution);
        }
    }
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    /* 
    choose possible heuristic function you want
    Manhattan, Euclidean, Diagonal, or 0 (Dijkstra)
    Remember tie_breaker learned in lecture, add it here ?
    *
    *
    *
    STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    please write your code below
    */
   //Dijkstra 没有贪心项
   double h = 0;
   //Manhattan 只走x,y,z轴   大于h*
    //h = abs(node2->coord(0) - node1->coord(0)) + abs(node2->coord(1) - node1->coord(1)) + abs(node2->coord(2) - node1->coord(2));

    //Eucilidean   两点之间最短距离  小于h*
    //h = sqrt(pow(node2->coord(0) - node1->coord(0),2) + pow(node2->coord(1) - node1->coord(1),2) + pow(node2->coord(2) - node1->coord(2),2));
    
    //diagonal 走对角线   应该是最佳h*
    
    double deltaX = abs(node1->coord(0) - node2->coord(0));
    double deltaY = abs(node1->coord(1) - node2->coord(1));
    double deltaZ = abs(node1->coord(2) - node2->coord(2));
    double minimum = min(deltaX, min(deltaY, deltaZ));
    double maximum = max(deltaX, max(deltaY, deltaZ));
    double middle = deltaX+deltaY+deltaZ-minimum-maximum;
    h = sqrt(3) * minimum + sqrt(2) * (middle-minimum) + maximum - middle;
    
    return h;

    //tie breaker: 对两条相同的路径，找到倾向性
    //return h*(1.0+1.0/100000);
}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt)
{   
    ros::Time time_1 = ros::Time::now();    

    //index of start_point and end_point  
    //起点与终点坐标
    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx   = coord2gridIndex(end_pt);
    goalIdx = end_idx;

    //position of start_point and end_point
    start_pt = gridIndex2coord(start_idx);
    end_pt   = gridIndex2coord(end_idx);

    //Initialize the pointers of struct GridNode which represent start node and goal node
    GridNodePtr startPtr = new GridNode(start_idx, start_pt);
    GridNodePtr endPtr   = new GridNode(end_idx,   end_pt);

    //openSet is the open_list implemented through multimap in STL library
    openSet.clear();
    // currentPtr represents the node with lowest f(n) in the open_list
    GridNodePtr currentPtr  = NULL;
    GridNodePtr neighborPtr = NULL;

    //put start node in open set  初始点放进open set
    startPtr -> gScore = 0;
    startPtr -> fScore = getHeu(startPtr,endPtr);   
    //STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    startPtr -> id = 1; 
    startPtr -> coord = start_pt;
    //multimap的类型是<double, GridNode*> 这里默认是升序排列
    openSet.insert( make_pair(startPtr -> fScore, startPtr) );
    /*
    STEP 2 :  some else preparatory works which should be done before while loop
    please write your code below
    */
    
    vector<GridNodePtr> neighborPtrSets;
    vector<double> edgeCostSets;

    // this is the main loop
    while ( !openSet.empty() ){
        /*
        *
        *
        step 3: Remove the node with lowest cost function from open set to closed set
        please write your code below
        
        IMPORTANT NOTE!!!
        This part you should use the C++ STL: multimap, more details can be find in Homework description
        *
        *
        */
       //因为是升序，取表头的最小值，然后从openlist里放到closedlist里
        multimap<double, GridNodePtr> ::iterator it = openSet.begin();
        Vector3i cur_index = it->second->index;
        currentPtr = GridNodeMap[cur_index(0)][cur_index(1)][cur_index(2)];
        currentPtr->id = -1;
        openSet.erase(it);
        
        // if the current node is the goal  
        //找到目标点，结束
        if( currentPtr->index == goalIdx ){
            ros::Time time_2 = ros::Time::now();
            terminatePtr = currentPtr;
            ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );            
            return;
        }
        //get the succetion  只是获取邻点和边的代价
        AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);  //STEP 4: finish AstarPathFinder::AstarGetSucc yourself     

        /*
        *
        *
        STEP 5:  For all unexpanded neigbors "m" of node "n", please finish this for loop
        please write your code below
        *        
        */         
        for(int i = 0; i < (int)neighborPtrSets.size(); i++){
            /*
            *
            *
            Judge if the neigbors have been expanded
            please write your code below
            
            IMPORTANT NOTE!!!
            neighborPtrSets[i]->id = -1 : expanded, equal to this node is in close set
            neighborPtrSets[i]->id = 1 : unexpanded, equal to this node is in open set
            neighborPtrSets[i]->id = 0 : unexpanded, 既不在open set 也不在 closed set
            *        
            */
            neighborPtr = neighborPtrSets[i];

            if(neighborPtr -> id == 0){ //discover a new node, which is not in the closed set and open set
                /*
                *
                *
                STEP 6:  As for a new node, do what you need do ,and then put neighbor in open set and record it
                please write your code below
                *        
                */
               //更新g和f
                neighborPtr -> gScore = currentPtr -> gScore + edgeCostSets[i];
                neighborPtr -> fScore  = neighborPtr -> gScore + getHeu(neighborPtr, endPtr);
                //加入open list
                neighborPtr->id = 1;
                neighborPtr->cameFrom = currentPtr;
                openSet.insert(std::make_pair(neighborPtr -> fScore, neighborPtr));
                //continue;
            }
            else if(neighborPtr -> id == 1){ //this node is in open set and need to judge if it needs to update, the "0" should be deleted when you are coding
                /*
                *
                *
                STEP 7:  As for a node in open set, update it , maintain the openset ,and then put neighbor in open set and record it
                please write your code below
                *        
                */
               //更新g和f
                double new_g = currentPtr -> gScore + edgeCostSets[i];
                if(new_g < neighborPtr->gScore) 
                {
                    neighborPtr->gScore = new_g; 
                    neighborPtr -> fScore  = neighborPtr -> gScore + getHeu(neighborPtr, endPtr);
                    //更新父节点
                     neighborPtr -> cameFrom = currentPtr;
                     //还要在multimap中更新key值，可是无法把之前那个key删除？openSet会越来越大？
                     openSet.insert(std::make_pair(neighborPtr -> fScore, neighborPtr));
                }
                //continue;
            }
            else{//this node is in closed set
                /*
                *
                please write your code below
                *        不用处理
                */
                continue;
            }
        }      
    }
    //if search fails
    ros::Time time_2 = ros::Time::now();
    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
}


vector<Vector3d> AstarPathFinder::getPath() 
{   
    vector<Vector3d> path;
    vector<GridNodePtr> gridPath;
    /*
    *
    *
    STEP 8:  trace back from the curretnt nodePtr to get all nodes along the path
    please write your code below
    *      
    */
   GridNodePtr cur = terminatePtr;
   while(cur!=nullptr)
   {
       gridPath.push_back(GridNodeMap[cur->index(0)][cur->index(1)][cur->index(2)]);
       cur = cur->cameFrom;
   }
    for (auto ptr: gridPath)
        path.push_back(ptr->coord);
        
    reverse(path.begin(),path.end());

    return path;
}