import time 
# from tqdm import tqdm
from progress.bar import Bar

time.sleep(2)
# pbar = tqdm(total=100)
bar = Bar('Processing', max=20)
for i in range(20):
    # Do some work
    bar.next()
bar.finish()

time.sleep(2)
# for i in range(100):
#         pbar.update(1)
#         time.sleep(0.1)
#         # pbar.set_postfix(
#         #     epoch=batch_idx,
#         #     idx_train=len(test_loader)*BATCH_SIZE,
#         #     loss=f"{batch_loss:.6f}",
#         #     accuracy=f"{batch_acc:.6f}",
#         #     J_Precision=f"{J_TP/(J_TP+J_FP+1.e-10):.6f}",
#         #     J_Recall=f"{J_TP/(J_TP+J_FN+1.e-10):.6f}"
#         # )


# pbar.close()

