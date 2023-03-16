# COVID-19 model for Sweden


To run the model use the command

```matlab
solve_Swedish_outbreak([1:5],'ResultFileBaseName','')
```
### Parameters
- The first parameter specifies which of the contact rate scenarios defined in `solve_Swedish_outbreak.m` are run.
- Second parameter sets the base filename for each resulting `.mat`-file.
- The third parameter can take an optional struct which let's you configure specific parameter values for the model. Refer to `./model/parameters_Swe_Corona_Radiation.m` for what these are.