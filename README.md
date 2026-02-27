# Overdose Prevention Centers and Neighborhood Commercial Activity in New York City
*Bennett Allen, PhD; Cale Basaraba, MPH; Laura C. Chambers, PhD; Czarina N. Behrends, PhD; Brandon D.L. Marshall, PhD; Magdalena Cerda, DrPH*

This is a code repository for the above manuscript:

https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2845604

Script descriptions:
- `01_make_isochrones.R` creates 5- and 10-minute walking buffers around all treated and donor sites.
- `02_aggregate_buffer_data.R` aggregates SafeGraph foot traffic and consumer spending data within each site buffer.
- `03_biweekly_covariate_setup.R` sets up SafeGraph analysis at a biweekly timescale and integrates covariates.
- `04_safegraph_models.R` runs all augmented synthetic control models presented in the manuscript.

No data is included in this repository.

All code inquiries: cale.basaraba@nyulangone.org
