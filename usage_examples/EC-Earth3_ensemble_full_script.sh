#!/usr/bin/env bash
###IPSL had 15 historical members downloaded and 10 future, so this is the pipeline for that.
### We're assuming you're in the right dir already to run the core_scripts


model="EC-Earth3"

for i in {101..127}; do
   python precip_projection.py --model $model --member r${i}i1p1f1 --experiment historical
done


for i in {101..127}; do
    python precursor_projection.py --model $model --member r${i}i1p1f1 --experiment historical --variables z500 u850 v850
done

for i in {101..109} {115..120}; do
    python precip_projection.py --model $model --member r${i}i1p1f1 --experiment ssp370
    python z500_detrend.py --model $model --member r${i}i1p1f1 --experiment ssp370
    python precursor_projection.py --model $model --member r${i}i1p1f1 --experiment ssp370 --variables z500_detrend u850 v850

done

chmod -R g+rw /Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/raw/${model}/z500_detrend/ssp370/
chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/indices/${model}/

# decomposition computed using full ensemble estimate
bash run_decomposition.sh --model $model \
--members "r101i1p1f1 r102i1p1f1 r103i1p1f1 r104i1p1f1 r105i1p1f1 r106i1p1f1 r107i1p1f1 r108i1p1f1 \
r109i1p1f1 r110i1p1f1 r111i1p1f1 r112i1p1f1 r113i1p1f1 r114i1p1f1 r115i1p1f1 r116i1p1f1 r117i1p1f1 \
r118i1p1f1 r119i1p1f1 r120i1p1f1 r121i1p1f1 r122i1p1f1 r123i1p1f1 r124i1p1f1 r125i1p1f1 r126i1p1f1 r127i1p1f1" \
--future_members "r101i1p1f1 r102i1p1f1 r103i1p1f1 r104i1p1f1 r105i1p1f1 r106i1p1f1 r107i1p1f1 \
r108i1p1f1 r109i1p1f1 r115i1p1f1 r116i1p1f1 r117i1p1f1 r118i1p1f1 r119i1p1f1 r120i1p1f1" \
--overwrite


#individual member-to-member decomposition:
for i in {101..109} {115..120}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f1 
done

#a bit of duplication here to get biases for the last 5 hist members:
for i in {101..127}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f1 \
    --futureexp "none"
done

SAVEDIR="/Data/skd/projects/global/cmip6_precursors/outputs/"

#Now we compute trend terms for each future member based on the historical ensemble, and we don't need
#to run the run_decomposition.sh script as the decomps are already stored now.
for i in {101..109} {115..120}; do
    #uncertainty estimate for trend based on ens mean decomp:
    python compute_terms.py \
    --ref_model ERA5 --hist_model $model --future_model $model \
    --hist_experiment historical --future_experiment ssp370 \
    --hist_member ens --future_member r${i}i1p1f1 \
    --savedir $SAVEDIR --overwrite
done

chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/decompositions/${model}/
chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/terms/${model}/
