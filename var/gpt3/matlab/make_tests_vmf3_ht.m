function [] = make_tests_vmf3_ht (num_tests)

% load the statistics package (for unifrnd)
pkg load statistics

printf("[[maybe_unused]] constexpr const int num_tests = %d;\n", num_tests);

mjd = unifrnd(44239e0, 62502e0, num_tests, 1); % random MJD in range 1980-2030
% printf("dso::modified_julian_day mjd{%d};\n", floor(mjd));

% mjd fraction to nanoseconds
printf("const dso::datetime<dso::nanoseconds> ts[] = {\n");
for i=1:num_tests-1
	frac = mjd(i) - floor(mjd(i));
	%printf("dso::nanoseconds nanosec{%dL};\n", int64(frac*86400e0*1e9));
	% c++ datetime
	%printf("dso::datetime<dso::nanoseconds> t{mjd, nanosec};\n");
	printf("\tdso::datetime<dso::nanoseconds>(dso::modified_julian_day(%d), dso::nanoseconds(%dL)),\n", floor(mjd(i)), int64(frac*86400e0*1e9));
end
frac = mjd(num_tests) - floor(mjd(num_tests));
printf("\tdso::datetime<dso::nanoseconds>(dso::modified_julian_day(%d), dso::nanoseconds(%dL))};\n", floor(mjd(num_tests)), int64(frac*86400e0*1e9));

lat = unifrnd(-pi/2, pi/2, num_tests, 1); % random latitude vector
printf("const double lats[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", lat(i));
end
printf("%+.20f};\n", lat(num_tests));

lon = unifrnd(-pi, pi, num_tests, 1); % random longitude vector
printf("const double lons[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", lon(i));
end
printf("%+.20f};\n", lon(num_tests));

h_ell = unifrnd(-100e0, 3e3, num_tests, 1); %  random ellipsoidal height vector
printf("const double hgts[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", h_ell(i));
end
printf("%+.20f};\n", h_ell(num_tests));

zd = unifrnd(1e-1, pi/2, num_tests, 1); %  random zenith angles
printf("const double zds[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", zd(i));
end
printf("%+.20f};\n", zd(num_tests));

ah = unifrnd(1e-3, 2e-2, num_tests, 1); %  random ah values
printf("const double ahs[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", ah(i));
end
printf("%+.20f};\n", ah(num_tests));

aw = unifrnd(1e-4, 2e-3, num_tests, 1); %  random aw values
printf("const double aws[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", aw(i));
end
printf("%+.20f};\n", aw(num_tests));

mfhr = zeros(num_tests, 1);
mfwr = zeros(num_tests, 1);
for i=1:num_tests
	[ mfhr(i) , mfwr(i) ] = vmf3_ht ( ah(i) , aw(i) , mjd(i) , lat(i) , lon(i) , h_ell(i) , zd(i) );
end

printf("const double mfhr[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", mfhr(i));
end
printf("%+.20f};\n", mfhr(num_tests));

printf("const double mfwr[] = {\n");
for i=1:num_tests-1
	printf("\t%+.20f,\n", mfwr(i));
end
printf("%+.20f};\n", mfwr(num_tests));

printf("double mfh,mfw;\n");
printf("for (int i=0; i<%d; i++) {\n", num_tests);
%for i=1:num_tests-1
%    j = i - 1;
printf("\tdso::vmf3_ht(ahs[i], aws[i], ts[i], lats[i], lons[i], hgts[i], zds[i], mfh, mfw);\n");
printf("\tprintf(\"DeltaMfh=%%+15.12e DeltaMfw=%%+15.12e\\n\", std::abs(mfh-mfhr[i]), std::abs(mfw-mfwr[i]));\n");
%end
printf("\n}\nreturn 0;\n}");
	
