import csv
import pandas as pd
  
results = pd.read_csv('fundamentals.csv')
length = len(results)
 
# opening the CSV file
with open('fundamentals.csv', mode ='r')as file:
   
  # reading the CSV file
  csvFile = csv.reader(file)
  i = 0
  
  with open("nystock.txt", "w") as log:
      for lines in csvFile:
          if i == 0:
              i += 1
              continue
          
          data = lines[3:10]
          valid = True
          
          for d in data:
              if d == "":
                  valid = False
          
          if not valid:
              continue
          
          for d in data:
              if d == "":
                  log.write(str(0))
              else:
                  log.write(str(float(d)))
              log.write(" ")
          
          if i != length:
              log.write("\n")
          i+=1
          