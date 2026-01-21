library('Gviz')


true <- AnnotationTrack("groundtruth/GHIST_2025_single_sweep_final_goldstandard.bed", name="True Sweeps", fill="darkgreen")
pred <- AnnotationTrack("final_submissions/single_sweep/single_sweep_3446506.bed", name="Predicted Sweeps", fill="red")

plotTracks(list(true, pred),
           chromosome="15",
           from=3.5e7, to=3.6e7)

