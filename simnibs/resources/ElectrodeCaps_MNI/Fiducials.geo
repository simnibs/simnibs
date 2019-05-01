View "projected positions" {
SP(0.000000, 82.900000, -43.000000){2};
SP(0.000000, -116.200000, -30.500000){2};
SP(-79.600000, -20.600000, -48.600000){2};
SP(80.600000, -20.600000, -48.100000){2};
};
myView = PostProcessing.NbViews-1; // indexing starts with 0
View[myView].PointType=1; // spheres
View[myView].PointSize=7;

View "projected positions - labels" {
T3(0.000000, 82.900000, -43.000000, 0){"Nz"};
T3(0.000000, -116.200000, -30.500000, 0){"Iz"};
T3(-79.600000, -20.600000, -48.600000, 0){"LPA"};
T3(80.600000, -20.600000, -48.100000, 0){"RPA"};
};

