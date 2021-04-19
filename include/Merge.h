/*
 * @Description: 
 * @Author: CrazyMax
 * @Date: 2021-02-01 10:15:58
 * @LastEditTime: 2021-04-19 16:41:01
 * @LastEditors: CrazyMax
 */

#ifndef MERGE_H
#define MERGE_H

#include <list>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include "EDLines.h"
#include <time.h>
// multithreading
//#include <thread> 
//#include <mutex>
//#include <future>

using namespace std;
//std::mutex mtx;

typedef struct 
{
    double length;
    double angular;
    cv::Point2d start;
	cv::Point2d end;
    int sign = 0;
}LINEINFO;

class Merge 
{
public:
     /**
     * @description: fusion into real line segments
     * @param  {const vector<LS>} lines set
     * @param  {double} tau_theta angular proximity threshold
     * @param  {double} xi_s spatial proximity threshold
     * @return {vector<LS>}
     */
    vector<LS> MergeLine(const vector<LS> lines, double tau_theta, double xi_s); 

private:
    /**
     * @description: Get point to line foot point
     * @param  {const cv::Point2f} point out of line
     * @param {const cv::Point2f} line begin point
     * @param {const cv::Point2f} line end point
     * @return {double}
     */
    cv::Point2f GetLineDist(const cv::Point2d &pt, const cv::Point2d &begin, const cv::Point2d &end);

    /**
     * @description: compute line gradient
     * @param  {LS} line end-points
     * @return {vector<double>} format g dx dy
     */
    vector<double> LineGradient(LS endPoints);

    /**
     * @description: compute point distance
     * @param {Point2d} point 1
     * @param {Point2d} point 2
     * @return {double}
     */
    inline double PointDistance(cv::Point2d p1, cv::Point2d p2);

    /**
     * @description: 
     * @param {LS} endPoints input line end-points
     * @param {LINEINFO} &li return line length/angular/end-points
     * @return {*}
     */
    inline LINEINFO CalculationLineInfo(LS endPoints);

    /**
     * @description: denscending order quick sort line 
     * @param {vector<LINEINFO>} &lf_v lines set
     * @param {int} lo sort begin index
     * @param {int} hi sort end index
     */
    void QuickSort(vector<LINEINFO> &lf_v, int lo, int hi);
 
    /**
     * @description: merge two lines
     * @param {LINEINFO} L1 line 1
     * @param {LINEINFO} L2 line 2
     * @param {double} xi_s spatial proximity threshold
     * @param {double} tau_theta angular proximity threshold
     * @param {LINEINFO} &merge_li merge result
     */
    bool MergeTwoLine(LINEINFO L1, LINEINFO L2, double xi_s, double tau_theta, LINEINFO &merge_li);
};
#endif