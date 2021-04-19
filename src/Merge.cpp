/*
* @Description: 
* @Author: CrazyMax
* @Date: 2021-02-01 10:15:58
 * @LastEditTime: 2021-04-19 16:32:23
 * @LastEditors: CrazyMax
*/

#include <Merge.h>

cv::Point2f Merge::GetLineDist(const cv::Point2d &pt, const cv::Point2d &begin, const cv::Point2d &end)   
{
    cv::Point2f retVal;

    double dx = begin.x - end.x;
    double dy = begin.y - end.y;
    if(abs(dx) < 0.00000001 && abs(dy) < 0.00000001){
        retVal = begin;
        return retVal;
    }

    double u = (pt.x - begin.x)*(begin.x - end.x) +
                (pt.y - begin.y)*(begin.y - end.y);
    u = u/((dx*dx)+(dy*dy));

    retVal.x = begin.x + u*dx;
    retVal.y = begin.y + u*dy;

    return retVal;
}

vector<double> Merge::LineGradient(LS endPoints){
    //format: g dx dy
    vector<double> grad_d_v;
    double dx = endPoints.end.x - endPoints.start.x;
    double dy = endPoints.end.y - endPoints.start.y;
    if(dx){
        double g = dy / dx;
    }else{
        double g = pow(10, 10);
    }
    grad_d_v.push_back(dx);
    grad_d_v.push_back(dy);
    return grad_d_v;
}

inline double Merge::PointDistance(cv::Point2d p1, cv::Point2d p2){
    return sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));
}

inline LINEINFO Merge::CalculationLineInfo(LS endPoints){
    LINEINFO li;
    li.length = sqrt(pow((endPoints.end.x - endPoints.start.x), 2) + 
                    pow((endPoints.end.y - endPoints.start.y), 2));
    //https://blog.csdn.net/chelen_jak/article/details/80518973
    li.angular = atan2((endPoints.end.y - endPoints.start.y), (endPoints.end.x - endPoints.start.x)) ;             
    if(li.angular <= 0.0) li.angular += M_PI;//convert from [-pi, pi] to [0, pi]
    //li.angular = angle * (180 / M_PI);
    li.start = endPoints.start;
    li.end = endPoints.end;
    return li;
}

void Merge::QuickSort(vector<LINEINFO> &lf_v, int lo, int hi){
    int i = lo, j = hi;
    LINEINFO temp;
    if(i < j){
        temp = lf_v[i];
        while (i != j){
            while(j > i && lf_v[j].length <= temp.length)-- j;
            lf_v[i] = lf_v[j];
            while(i < j && lf_v[i].length >= temp.length)++ i;
            lf_v[j] = lf_v[i];
        }
        lf_v[i] = temp;
        QuickSort(lf_v, lo, i - 1);
        QuickSort(lf_v, i + 1, hi);
    }
}

bool Merge::MergeTwoLine(LINEINFO L1, LINEINFO L2, double xi_s, double tau_theta, LINEINFO &merge_li){
    vector<cv::Point2d> l1{L1.start, L1.end}, l2{L2.start, L2.end};
    // **** 1.if L2 length long than L1 swap than ****
    if(L2.length > L1.length)
    {
        swap(L1, L2);
        swap(l1, l2);
    }
    // **** 1.search min distance of two lines point ****
    double min_dis[2][2];
    int min_i = 0, min_j =0;
    double min_norm_value = norm(l1[0] - l2[0]);
    //use multithreading compute
    //future<void> ft1 = async(std::launch::async, [&]{
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                //unique_lock<mutex> lock(mtx);
                min_dis[i][j] = norm(l1[i] - l2[j]);
                if(min_dis[i][j] < min_norm_value){
                    min_i = i;
                    min_j = j;
                    min_norm_value = min_dis[i][j];
                }
            }
        }   
    //});
        
    // **** 2.get pairwise distances between the 4 end-points ****
    vector<cv::Point2d> c, l; 
    c.push_back(l1[min_i]);
    c.push_back(l1[1 - min_i]);
    l.push_back(l2[min_j]);
    l.push_back(l2[1 - min_j]);

    double tau_s = xi_s * L1.length;  
    if(min_norm_value > tau_s){
        return false;
    }else{
        double lambda = (L2.length / L1.length) + (min_norm_value / tau_s);
        double tau_theta_star = (1.0 - (1.0 / ( 1.0 + exp(-2.0 * (lambda - 1.5))))) * tau_theta;
        if((fabs(L1.angular - L2.angular) < tau_theta_star) || (fabs(L1.angular - L2.angular) > (M_PI - tau_theta_star))){              
            //l1-1 & l1 - 2  l1-1 & l2 - 2   l1-2 & l2 - 1  l1-2 & l2 - 2
            double max_dis[2][2];
            double max_dist_value = 0.0;
            int max_i = 0, max_j =0;
                //use multithreading compute
                //future<void> ft2 = async(std::launch::async, [&]{
                for(int i = 0; i < 2; i++){
                    for(int j = 0; j < 2; j++){
                        //unique_lock<mutex> lock(mtx);
                        if(i == 0  && j == 0){
                            max_dis[i][j] = PointDistance(c[i], c[i + 1]);//PointDistance(c1, l1o);                          
                        }else{
                            max_dis[i][j] = PointDistance(c[i], l[j]);
                        }
                        if(max_dis[i][j] > max_dist_value){
                            max_i = i;
                            max_j = j;
                            max_dist_value = max_dis[i][j];
                        }
                    }
                }
            //});
            
            //gardient of merged line shold be similar to gradient of line L1
            //final check on the absolute angular difference of the longer segment L1 and the merged segment M               
            //final check on the absolute angular difference of the longer segment L1 and the merged segment M.
            if(max_i == 0 && max_j == 0){
                LS merge_ls(c[max_i], c[1- max_i]);
                merge_li = CalculationLineInfo(merge_ls);
            }else{
                LS merge_ls(c[max_i], l[max_j]);
                merge_li = CalculationLineInfo(merge_ls);
            }
            if(fabs(L1.angular - merge_li.angular) < (tau_theta / 2)) return true;
            else return false;
        }
        else{
            return false;
        }           
    }           
}

vector<LS> Merge::MergeLine(const vector<LS> lines, double tau_theta, double xi_s){
    unordered_multimap<int, LINEINFO> uMap;
    vector<LINEINFO> lf_v;
    for(auto line:lines)
        lf_v.push_back(CalculationLineInfo(line));
    // **** record original line number ****
    size_t num_lines = lf_v.size();
    while(1){
        // **** 1.sort line abort line length ****
        size_t lf_v_size = lf_v.size();
        // sort with line length
        QuickSort(lf_v, 0, (int)lf_v_size - 1);
        uMap.clear();
        for(int index = 0; index < lf_v.size(); index++)
            uMap.insert({index, lf_v[index]});
        int i = 0;
        for(auto it = uMap.begin(); it != uMap.end(); ++it){
            // **** 2.find lines with direction close to L1 ****
                vector<int> angular_sign_num;
                for(auto it_a = it; it_a != uMap.end(); ++it_a){
                    if(it_a->second.length == it->second.length) continue;
                    if(fabs(it_a->second.angular - it->second.angular) < tau_theta){
                        angular_sign_num.push_back(it_a->first);
                    }
                }
                if(angular_sign_num.size() > 0){
                        double tau_s = xi_s * it->second.length;
                        //TODO -- while line between middle but near to long line, can no to merge, because distance long than set
                        //TODO -- some short line in corner can no to merge
                        vector<int> spatial_sign_num;
                        for(int i_x = 0; i_x < angular_sign_num.size(); i_x++){
                            //judge L1 x diff with other start and end x
                            //it->second == L1
                            unordered_multimap<int, LINEINFO>::iterator iter_sp = uMap.find(angular_sign_num[i_x]);
                            if(iter_sp != uMap.end()){
                                if(fabs(it->second.start.x - iter_sp->second.start.x) < tau_s || fabs(it->second.start.x - iter_sp->second.end.x) < tau_s ||
                                fabs(it->second.end.x - iter_sp->second.start.x) < tau_s || fabs(it->second.end.x - iter_sp->second.end.x) < tau_s){
                                //judge L1 y diff with other start and end y
                                if(fabs(it->second.start.y - iter_sp->second.start.y) < tau_s || fabs(it->second.start.y - iter_sp->second.end.y) < tau_s ||
                                    fabs(it->second.end.y - iter_sp->second.start.y) < tau_s || fabs(it->second.end.y - iter_sp->second.end.y) < tau_s){
                                    spatial_sign_num.push_back(iter_sp->first);
                                    }
                                }
                            }          
                        }

                        // **** 4.merge short line L2 to long line L1 ****
                        for(int i_m = 0; i_m < spatial_sign_num.size(); i_m++){
                            //iter_mg->second == L2
                            unordered_multimap<int, LINEINFO>::iterator iter_mg = uMap.find(spatial_sign_num[i_m]);
                            if(iter_mg != uMap.end()){                        
                                //mergeTowLine
                                LINEINFO merge_li;
                                // TODO merge two line use a lot time
                                bool is_merge = MergeTwoLine(it->second, iter_mg->second, xi_s, tau_theta, merge_li);
                                if(is_merge){
                                    //use merge line replace original line
                                    it->second = merge_li;    
                                    iter_mg->second.length = 9999.0;// mark bad line           
                                }                           
                            }
                        }
                        // **** 5.remove have been merge to L1 lines **** 
                        for(auto it_up = uMap.begin(); it_up != uMap.end(); ++it_up){
                            if(it_up->second.length == 9999.0){
                            uMap.erase(it_up->first);                
                            }
                        }     
                    }
                }
                lf_v.clear();
                for(auto it_up = uMap.begin(); it_up != uMap.end(); ++it_up){
                    lf_v.push_back(it_up->second);
            }
        //judge if have any line have been merge
        if(lf_v.size() == num_lines) break;
        else num_lines = lf_v.size();
    }
    vector<LS> t_ls_v;
    for(auto lf_info:lf_v){
        t_ls_v.push_back(LS(lf_info.start, lf_info.end));
    }
    return t_ls_v;
}  