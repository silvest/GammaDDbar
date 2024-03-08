/* 
 * File:   histo.h
 * Author: silvest
 *
 * Created on January 21, 2014, 12:00 PM
 */

#ifndef HISTO_H
#define	HISTO_H
#include <string>
#include <vector>
#include <map>
#include <TH1D.h>
#include <TH2D.h>

using namespace std;

class histo {
public:
    histo(const histo& orig);
    virtual ~histo() {};
    map<string, TH1D*> h1d;
    map<string, TH2D*> h2d;
    vector<string> h1dnames;
    vector<string> h2dnames;

    histo(map<string, double>& obs);

    void createH1D(string name, int binx, double minx, double maxx);

    void fillh1d();

    void createH2D(string namey, string namex, int binx, double minx, double maxx,
            int biny, double miny, double maxy);

    void fillh2d();

    void write();
    map<string, double>& myobs;
private:

};

#endif	/* HISTO_H */

