function [] = make_tests_31 (num_tests)

% load the statistics package (for unifrnd)
pkg load statistics

% load the grid file
grid = gpt3_1_fast_readGrid();

mjd = unifrnd(44239e0, 62502e0); % random MJD in range 1980-2030
printf("dso::modified_julian_day mjd{%d};\n", floor(mjd));

% mjd fraction to nanoseconds
frac = mjd - floor(mjd);
printf("dso::nanoseconds nanosec{%dL};\n", int64(frac*86400e0*1e9));

% c++ datetime
printf("dso::datetime<dso::nanoseconds> t{mjd, nanosec};\n");

lat = unifrnd(-pi/2, pi/2, num_tests); % random latitude vector
printf("const double lats[] = {");
for i=1:num_tests-1
	printf("%+.20f,", lat(i));
end
printf("%+.20f};\n", lat(num_tests));

lon = unifrnd(-pi, pi, num_tests); % random longitude vector
printf("const double lons[] = {");
for i=1:num_tests-1
	printf("%+.20f,", lon(i));
end
printf("%+.20f};\n", lon(num_tests));

h_ell = unifrnd(-100e0, 3e3, num_tests); %  random ellipsoidal height vector
printf("const double hgts[] = {");
for i=1:num_tests-1
	printf("%+.20f,", h_ell(i));
end
printf("%+.20f};\n", h_ell(num_tests));

[p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w] = gpt3_1_fast (mjd,lat,lon,h_ell,0,grid);

printf("dso::gpt3_result gout[%d];\n", num_tests);
printf("const dso::gpt3_result res35it0[] = {");
for i=1:num_tests-1
	printf("{%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f},\n", 
		p(i),T(i),dT(i),Tm(i),e(i),ah(i),aw(i),la(i),undu(i),Gn_h(i),Ge_h(i),Gn_w(i),Ge_w(i));
end
printf("{%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f}};\n", 
        p(num_tests),T(num_tests),dT(num_tests),Tm(num_tests),e(num_tests),ah(num_tests),aw(num_tests),la(num_tests),undu(num_tests),Gn_h(num_tests),Ge_h(num_tests),Gn_w(num_tests),Ge_w(num_tests));

printf("if (dso::gpt3_fast(t, lats, lons, hgts, %d, 0, &grid35, gout)) {\n", num_tests);
printf("  fprintf(stderr,""ERROR. Failed call to gpt3_fast\\n"");\n  return 5;\n");
printf("}\n");
printf("for (int i=0; i<%d; i++) {\n", num_tests);
printf("  assert(std::abs(gout[i].p-res35it0[i].p) < PRECISION);\n");
printf("  assert(std::abs(gout[i].T-res35it0[i].T) < PRECISION);\n");
printf("  assert(std::abs(gout[i].dT-res35it0[i].dT) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Tm-res35it0[i].Tm) < PRECISION);\n");
printf("  assert(std::abs(gout[i].e-res35it0[i].e) < PRECISION);\n");
printf("  assert(std::abs(gout[i].ah-res35it0[i].ah) < PRECISION);\n");
printf("  assert(std::abs(gout[i].aw-res35it0[i].aw) < PRECISION);\n");
printf("  assert(std::abs(gout[i].la-res35it0[i].la) < PRECISION);\n");
printf("  assert(std::abs(gout[i].undu-res35it0[i].undu) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Gn_h-res35it0[i].Gn_h) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Ge_h-res35it0[i].Ge_h) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Gn_w-res35it0[i].Gn_w) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Ge_w-res35it0[i].Ge_w) < PRECISION);\n");
printf("}\n");

[p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w] = gpt3_1_fast (mjd,lat,lon,h_ell,1,grid);

printf("const dso::gpt3_result res35it1[] = {");
for i=1:num_tests-1
	printf("{%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f},\n", 
		p(i),T(i),dT(i),Tm(i),e(i),ah(i),aw(i),la(i),undu(i),Gn_h(i),Ge_h(i),Gn_w(i),Ge_w(i));
end
printf("{%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f}};\n", 
        p(num_tests),T(num_tests),dT(num_tests),Tm(num_tests),e(num_tests),ah(num_tests),aw(num_tests),la(num_tests),undu(num_tests),Gn_h(num_tests),Ge_h(num_tests),Gn_w(num_tests),Ge_w(num_tests));

printf("if (dso::gpt3_fast(t, lats, lons, hgts, %d, 1, &grid35, gout)) {\n", num_tests);
printf("  fprintf(stderr,""ERROR. Failed call to gpt3_fast\\n"");\n  return 5;\n");
printf("}\n");
printf("for (int i=0; i<%d; i++) {\n", num_tests);
printf("  assert(std::abs(gout[i].p-res35it1[i].p) < PRECISION);\n");
printf("  assert(std::abs(gout[i].T-res35it1[i].T) < PRECISION);\n");
printf("  assert(std::abs(gout[i].dT-res35it1[i].dT) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Tm-res35it1[i].Tm) < PRECISION);\n");
printf("  assert(std::abs(gout[i].e-res35it1[i].e) < PRECISION);\n");
printf("  assert(std::abs(gout[i].ah-res35it1[i].ah) < PRECISION);\n");
printf("  assert(std::abs(gout[i].aw-res35it1[i].aw) < PRECISION);\n");
printf("  assert(std::abs(gout[i].la-res35it1[i].la) < PRECISION);\n");
printf("  assert(std::abs(gout[i].undu-res35it1[i].undu) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Gn_h-res35it1[i].Gn_h) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Ge_h-res35it1[i].Ge_h) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Gn_w-res35it1[i].Gn_w) < PRECISION);\n");
printf("  assert(std::abs(gout[i].Ge_w-res35it1[i].Ge_w) < PRECISION);\n");
printf("}\n");
