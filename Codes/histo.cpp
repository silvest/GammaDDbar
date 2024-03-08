#include "histo.h"
#include <iostream>

histo::histo(map<string, double>& obs) : h1d(), h2d(), myobs(obs) {
};

void histo::createH1D(string name, int binx, double minx, double maxx) {
//    std::cout << "creating h1d " << name << std::endl;
    h1dnames.push_back(name);
    h1d[name] = new TH1D(name.c_str(), name.c_str(), binx, minx, maxx);
}

void histo::fillh1d() {
    for (vector<string>::iterator it = h1dnames.begin(); it != h1dnames.end(); ++it) {
        string name = *it;
//        std::cout << "filling h1d " << name << std::endl;
        h1d[name]->Fill(myobs.at(name));
//        std::cout << "...done" << std::endl;
    }
}

void histo::createH2D(string namey, string namex, int binx, double minx, double maxx,
        int biny, double miny, double maxy) {
    string name = namey + "_vs_" + namex;
//    std::cout << "creating h2d " << name << std::endl;
    h2dnames.push_back(name);
    h2d[name] = new TH2D(name.c_str(), name.c_str(), binx, minx, maxx, biny, miny, maxy);
}

void histo::fillh2d() {
    for (vector<string>::iterator it = h2dnames.begin(); it != h2dnames.end(); ++it) {
        string name = *it;
//        std::cout << "filling h2d " << name << std::endl;
        unsigned pos = name.find("_vs_");
        double obsx = myobs.at(name.substr(pos + 4));
        double obsy = myobs.at(name.substr(0, pos));
        h2d.at(name)->Fill(obsx, obsy);
//        std::cout << "...done" << std::endl;
    }
}

void histo::write() {
        for (vector<string>::iterator it = h1dnames.begin(); it != h1dnames.end(); ++it) {
        string name = *it;
        h1d[name]->Write();
    }
        for (vector<string>::iterator it = h2dnames.begin(); it != h2dnames.end(); ++it) {
        string name = *it;
        h2d[name]->Write();
    }
}
