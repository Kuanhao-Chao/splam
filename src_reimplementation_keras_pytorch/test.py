# Example of target with class indices
import torch
import torch.nn as nn
loss = nn.CrossEntropyLoss()
input = torch.randn(3, 5, requires_grad=True)
print("input: ", input)
print("input: ", input.size())
target = torch.empty(3, dtype=torch.long).random_(5)
print("target: ", target)
print("target: ", target.size())

output = loss(input, target)
output.backward()
input = torch.randn(3, 5, requires_grad=True)

print("input: ", input)
print("input: ", input.size())
target = torch.randn(3, 5).softmax(dim=1)
output = loss(input, target)

print("output: ", output)
print("output: ", output.size())
output.backward()

