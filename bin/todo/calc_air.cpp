double calc_air (int active, vector<Residue> &res, vector<Coor3f *> &actpass) {
  const double upperlimit = 2;
  const double wellsize = 1;
  const double deffsum_lim = pow(upperlimit, -6);

  double deffsum = 0;
  vector< vector<double> > distances;
  Residue &actres = res[active];
  int n;
  for (n = 0; n < actres.atoms.size(); n++) {
    Coor3f &actatom = actres.atoms[n];
    vector<double> currdistances;
    //printf("n %d %d\n", n, active);
    for (int nn = 0; nn < actpass.size(); nn++) {
      //printf("nn %d\n", nn);
      Coor3f &passatom = *actpass[nn];
      double dissq = (actatom.x-passatom.x)*(actatom.x-passatom.x) +
       (actatom.y-passatom.y)*(actatom.y-passatom.y) +
       (actatom.z-passatom.z)*(actatom.z-passatom.z);
      double deffsum_contrib = pow(dissq * dissq * dissq, -1);
      //printf("contrib %.3e\n", deffsum_contrib);
      deffsum += deffsum_contrib;
      if (deffsum >= deffsum_lim) return 0;
      currdistances.push_back(dissq);
    }
    distances.push_back(currdistances);
  }

  double deff = pow(deffsum, -1.0/6);
  //printf("deff %.3f\n", deff);
  double force = force_max;
  if (deff < upperlimit + wellsize) force *= (deff - upperlimit);

  double gradfac = -1.0/6 * pow(deffsum, -7.0/6) * force;
  //printf("gradfac %.3e\n", gradfac);
  for (n = 0; n < actres.atoms.size(); n++) {
    Coor3f &actatom = actres.atoms[n];
    vector<double> &currdistances = distances[n];
    set<Coor3f *>::iterator it;
    int nn = 0;
    for (int nn = 0; nn < actpass.size(); nn++) {
      Coor3f &passatom = *actpass[nn];
      double &dissq = currdistances[nn];
      double disfac = 6 * pow(dissq, -4) * gradfac;
      double gradx = disfac * (actatom.x-passatom.x);
      double grady = disfac * (actatom.y-passatom.y);
      double gradz = disfac * (actatom.z-passatom.z);
      actatom.gradx += gradx; actatom.grady += grady; actatom.gradz += gradz;
      passatom.gradx -= gradx; passatom.grady -= grady; passatom.gradz -= gradz;
    }
  }
  double energy = 0;
  deff -= upperlimit;
  if (deff < wellsize) {
    energy = .5 * force_max * deff * deff;
  }
  else {
    energy = .5 * force_max * wellsize * wellsize + force_max * (deff - wellsize);
  }
  return energy;
}
