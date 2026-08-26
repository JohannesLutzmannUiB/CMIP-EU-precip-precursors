
###IPSL had 15 historical members downloaded and 10 future, so this is the pipeline for that.
### We're assuming you're in the right dir already to run the core_scripts


model="UKESM1-0-LL"

for i in {1..4} {8..10}; do
   python precip_projection.py --model $model --member r${i}i1p1f2 --experiment historical
done

for i in {5..7}; do
   python precip_projection.py --model $model --member r${i}i1p1f3 --experiment historical
done


for i in {1..4} {8..10}; do
    python precursor_projection.py --model $model --member r${i}i1p1f2 --experiment historical --variables z500 u850 v850
done

for i in {5..7}; do
    python precursor_projection.py --model $model --member r${i}i1p1f3 --experiment historical --variables z500 u850 v850
done

for i in {1..10}; do
    python precip_projection.py --model $model --member r${i}i1p1f2 --experiment ssp370
    python z500_detrend.py --model $model --member r${i}i1p1f2 --experiment ssp370
    python precursor_projection.py --model $model --member r${i}i1p1f2 --experiment ssp370 --variables z500_detrend u850 v850

done

chmod -R g+rw /Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/raw/${model}/z500_detrend/ssp370/
chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/indices/${model}/

# decomposition computed using full ensemble estimate
bash run_decomposition.sh --model $model \
--members "r1i1p1f2 r2i1p1f2 r3i1p1f2 r4i1p1f2 r5i1p1f3 r6i1p1f3 r7i1p1f3 r8i1p1f2 r9i1p1f2 r10i1p1f2" \
--future_members "r1i1p1f2 r2i1p1f2 r3i1p1f2 r4i1p1f2 r5i1p1f2 r6i1p1f2 r7i1p1f2 r8i1p1f2 r9i1p1f2 r10i1p1f2" \
--overwrite


#individual member-to-member decomposition:
for i in {1..4} {8..10}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f2 
done

#individual member-to-member decomposition:
for i in {5..7}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f3 \
    --future_members r${i}i1p1f2
done

#a bit of duplication here to get biases for the last 5 hist members:
for i in {1..4} {8..10}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f2 \
    --futureexp "none"
done


for i in {5..7}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f3 \
    --futureexp "none"
done

SAVEDIR="/Data/skd/projects/global/cmip6_precursors/outputs/"

#Now we compute trend terms for each future member based on the historical ensemble, and we don't need
#to run the run_decomposition.sh script as the decomps are already stored now.
for i in {1..10}; do
    #uncertainty estimate for trend based on ens mean decomp:
    python compute_terms.py \
    --ref_model ERA5 --hist_model $model --future_model $model \
    --hist_experiment historical --future_experiment ssp370 \
    --hist_member ens --future_member r${i}i1p1f2 \
    --savedir $SAVEDIR --overwrite
done

chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/decompositions/${model}/
chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/terms/${model}/
