import pandas as pd
import matplotlib.pyplot as plt

annotations = ["chess/", "gencode/", "refseq_ucsc/"]

for library in ["polyA", "ribozero"]:
    for annotation in annotations:
        before_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "BEFORE.tsv", delimiter="\t", index_col=0)
        after_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "AFTER.tsv", delimiter="\t", index_col=0)

        # before_df[]
        diff_df = after_df.iloc[[4]] - before_df.iloc[[4]] 
        print(diff_df)


        y = diff_df.values.tolist()[0]
        y = [float(value) for value in y]
        print("y: ", y)
        
        # Create a list of colors based on the values
        colors = ['blue' if value > 0 else 'red' for value in y]
        print("colors: ", colors)
        # Create the bar chart
        plt.bar(list(before_df.columns.values)[:12] , y[:12], color = colors)

        # Add labels and title
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.title('Bar Chart with Blue and Red Bars')

        # Display the chart
        plt.show()
        # print("before_df: ", before_df.iloc[[4]])
        # print("after_df : ", after_df.iloc[[4]])