
# based on nemo tabplot.py

import sys
import matplotlib.pyplot as plt

def plot_data(file_name):
    # Function to plot data from a file
    with open(file_name, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()
        # Parse data from lines, skipping lines starting with '#', and convert to float
        data = [[float(val) for val in line.split()] for line in lines if not line.startswith('#')]
        # Plot data using matplotlib
        for i in range(1, len(data[0])):
            plt.plot([row[0] for row in data], [row[i] for row in data])
        # Show the plot
        plt.show()

def main():
    # Check if filename is provided as command-line argument
    if len(sys.argv) < 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)
    
    # Get filename from command-line argument
    file_name = sys.argv[1]
    # Plot data from the specified file
    plot_data(file_name)

if __name__ == '__main__':
    main()

#the following code could be implemented using pyqt5 in the future, but this is more of a proof of concept. 
