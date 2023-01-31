import matplotlib.pyplot as plt
import csv
import os

def plotDiffusion():
    file = open("../output/diffusion.csv", mode='r')
    
    reader = csv.reader(file)
    
    timestep = []
    var = []
    for row in reader:
        if row[0] != "timestep":
            timestep.append(float(row[0]))
            var.append(float(row[1]))
    
    plt.plot(timestep,var)
    plt.xlabel("Time step")
    plt.ylabel("Diffusion Var(t)")
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    
    plt.savefig("../output/plot/diffusion.png")
    plt.clf()

def plotRDF():
    for file in os.listdir("../output/rdf"):
        reader = csv.reader(open("../output/rdf/" + os.fsdecode(file)))
        distance = []
        densities = []
        for row in reader:
            if row[0] != "distance":
                distance.append(float(row[0]))
                densities.append(float(row[1]))
        plt.plot(distance,densities, label = file.replace("time_", "").replace(".csv",""))
    plt.xlabel("distance")
    plt.ylabel("densities")
    plt.legend()
    plt.ylim(bottom=0)
    plt.xlim(left=0)

    plt.savefig("../output/plot/rdf.png")
    plt.clf()

plotDiffusion()
plotRDF()