
###IPSL had 15 historical members downloaded and 10 future, so this is the pipeline for that.
### We're assuming you're in the right dir already to run the core_scripts


model="CNRM-CM6-1"

for i in {1..6} {11..15} {17..25} {27..28}; do
   python precip_projection.py --model $model --member r${i}i1p1f2 --experiment historical
done


for i in {1..6} {11..15} {17..25} {27..28}; do
    python precursor_projection.py --model $model --member r${i}i1p1f2 --experiment historical --variables z500 u850 v850
done

for i in {1..1}; do
    python precip_projection.py --model $model --member r${i}i1p1f2 --experiment ssp370
done
for i in {1..2} {4..6}; do
    python z500_detrend.py --model $model --member r${i}i1p1f2 --experiment ssp370
    python precursor_projection.py --model $model --member r${i}i1p1f2 --experiment ssp370 --variables z500_detrend u850 v850

done

chmod -R g+rw /Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/raw/${model}/z500_detrend/ssp370/
chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/indices/${model}/

# decomposition computed using full ensemble estimate
bash run_decomposition.sh --model $model \
--members "r1i1p1f2 r2i1p1f2 r3i1p1f2 r4i1p1f2 r5i1p1f2 r6i1p1f2 r11i1p1f2 r12i1p1f2 r13i1p1f2 \
r14i1p1f2 r15i1p1f2 r17i1p1f2 r18i1p1f2 r19i1p1f2 r20i1p1f2 r21i1p1f2 r22i1p1f2 r23i1p1f2 \
r24i1p1f2 r25i1p1f2 r27i1p1f2 r28i1p1f2" \
--future_members "r1i1p1f2" \
--overwrite


#individual member-to-member decomposition:
for i in {1..1}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f2 
done

#a bit of duplication here to get biases for the last 5 hist members:
for i in {1..6} {11..15} {17..25} {27..28}; do
    bash run_decomposition.sh --model $model \
    --members r${i}i1p1f2 \
    --futureexp "none"
done

SAVEDIR="/Data/skd/projects/global/cmip6_precursors/outputs/"

#Now we compute trend terms for each future member based on the historical ensemble, and we don't need
#to run the run_decomposition.sh script as the decomps are already stored now.
for i in {1..1}; do
    #uncertainty estimate for trend based on ens mean decomp:
    python compute_terms.py \
    --ref_model ERA5 --hist_model $model --future_model $model \
    --hist_experiment historical --future_experiment ssp370 \
    --hist_member ens --future_member r${i}i1p1f2 \
    --savedir $SAVEDIR --overwrite
done

chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/decompositions/${model}/
chmod -R g+rw /Data/skd/projects/global/cmip6_precursors/outputs/terms/${model}/
