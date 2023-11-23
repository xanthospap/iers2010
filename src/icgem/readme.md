# ICGEM Files

We are currently ok with parsing ICGEM files when the file describes a gravity field 
(i.e. `product_type gravity_field`).
ICGEM files include a format specifier that can be `“icgem1.0”` or `“icgem2.0”`. 
Both are (internally) handled, without the user ever needed to know/retrieve 
the actual file version.

## Testing

Tests are inlcuded in the [unit_test](test/unit_tests) directory. 
A couple of dedicated python scripts are available in [script](test/script/) dir, 
to assist validating results.

According to specifications (see [documentation](#ICGEM-format-2023)), the parsing 
of version 1 files should be ok, **but** they are yet to be tested. Actually, i have 
tried with the file `EIGEN-5C.gfc`, but it seems that is breaks the rules (i.e. 
the line `dot    2    1 -.337000000000e-11 0.160600000000e-10 0.0000e+00 0.0000e+00` 
is not valid according to [documentation](#ICGEM-format-2023). Hence, i have not 
checked this version.

## Non-gravity ICGEM files

No other product type is tested or considered!

## Issues

### Order of Stokes coefficients

According to [documentation](#ICGEM-format-2023), the maximum degree of the Stokes 
coefficients included in the file is recorded in the header. However, this does not 
hold for the maximum order. Hence, it is quite difficult to assure that the requested 
number of Stokes coefficients to parse (requested by the user) will the same with 
the actual order parsed.
The best thing to do, is to get the actual degree/order by probing the `StokesCoeffs` 
instance used to store the coefs.
```
  Icgem gfc('the gfc file');
  StokesCoeffs stokes;
  datetime<nanoseconds> t(...);

  if (gfc.parse_data(degree, order, t, stokes)) {
    fprintf(stderr, "ERROR parsing file\n");
    return 1;
  }

  max_parsed_degree = stokes.max_degree();
  max_parsed_order  = stokes.max_order();
  // or assert degree/order if they must satisfy some condition(s)
```
 
# References 

<a name="ICGEM-format-2023"></a>Christoph Förste, Franz Barthelmes, and E. Sinem Ince, 
["The ICGEM-format"](http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf),
*Version February 2023*
