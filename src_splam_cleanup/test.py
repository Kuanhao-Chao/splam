import re

input_string = """
        Base level:    42.5     |    36.8    |
        Exon level:    41.2     |    65.7    |
"""

# Extract all numbers using regular expressions
numbers = re.findall(r"\d+\.\d+", input_string)

# Print the extracted numbers
for number in numbers:
    print(number)
