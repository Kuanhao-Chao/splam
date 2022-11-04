import torch
from memory_transformer_xl import MemoryTransformerXL

model = MemoryTransformerXL(
    num_tokens = 20000,
    dim = 1024,
    heads = 8,
    depth = 8,
    seq_len = 100,
    mem_len = 256,            # short term memory (the memory from transformer-xl)
    lmem_len = 256,           # long term memory (memory attention network attending to short term memory and hidden activations)
    mem_write_iters = 2,      # number of iterations of attention for writing to memory
    memory_layers = [6,7,8],  # which layers to use memory, only the later layers are actually needed
    num_mem_kv = 128,         # number of memory key/values, from All-attention paper

)

x1 = torch.randint(0, 20000, (100, 100))
print(x1.size())
logits1, mem1 = model(x1)

print("logits1: ", logits1.size())
print("mem1: ", mem1)


x2 = torch.randint(0, 20000, (100, 100))
print(x2.size())
logits2, mem2 = model(x2, memories = mem1)

print("logits2: ", logits2.size())
print("mem2: ", mem2)

# and so on with carrying over memories...