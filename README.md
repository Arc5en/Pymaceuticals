# Pymaceuticals

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - Add your analysis here.
#  

# In[ ]:


# This is the annotated version of the code I used in the matplotlib module challenge, found in the README file. The code was heavily inspired from TA Drew's speedrun.
# Most comments here will also explain the code in depth.


# In[4]:


# Preparing to use modules/libraries/packages to analyse collected data.
# Testing if python can read csv files with matplotlib properly in combination with pandas.

# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
full_mouse_df = pd.merge(study_results, mouse_metadata, on="Mouse ID", how="left")

# Display the data table for preview
full_mouse_df


# In[6]:


# Checking the number of mice.
len(full_mouse_df["Mouse ID"].unique())


# In[8]:


# Starting to clean the data to accurately represent a population's individuals in one instance.

# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mouse_ids = full_mouse_df[full_mouse_df.duplicated(subset = ["Mouse ID", "Timepoint"])]["Mouse ID"].unique()
duplicate_mouse_ids


# In[9]:


# Optional: Get all the data for the duplicate mouse ID. 
duplicate_mouse_dataset = full_mouse_df[full_mouse_df["Mouse ID"] == 'g989']
duplicate_mouse_dataset


# In[14]:


# Creating the dataframe with no duplicate mouses based on their ID.

# Create a clean DataFrame by dropping the duplicate mouse by its ID.
no_dupes_df = full_mouse_df[full_mouse_df["Mouse ID"].isin(duplicate_mouse_ids) == False]
no_dupes_df


# In[15]:


# Checking the number of mice in the clean DataFrame.
len(no_dupes_df["Mouse ID"].unique())


# ## Summary Statistics

# In[23]:


# Making a table showing summary stats such as mean and variance based on drug treatment using pandas dataframe.

# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.
means = no_dupes_df.groupby("Drug Regimen").mean()["Tumor Volume (mm3)"]
medians = no_dupes_df.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
variances = no_dupes_df.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
standard_deviations = no_dupes_df.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
sems = no_dupes_df.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]

summary_stats = pd.DataFrame({
    "Mean Tumor Volume (mm3)": means,
    "Median Tumor Volume (mm3)": medians,
    "Tumor Volume Variances": variances,
    "Tumor Volume Standard Deviations": standard_deviations,
    "Tumor Volume Standard Errors": sems
})
summary_stats


# In[24]:


# Making a table showing summary stats such as mean and variance based on drug treatment using aggregate function in pandas.

# Generate a summary statistics table of mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line.
stats_aggregated = no_dupes_df.groupby("Drug Regimen").agg({"Tumor Volume (mm3)": ["mean", "median", "var", "std", "sem"]})
stats_aggregated


# ## Bar and Pie Charts

# In[26]:


# Creating bar plot using pandas dataframes. Data used in this and subsequent plots are based on table created prior.

# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.

tests = no_dupes_df["Drug Regimen"].value_counts()
tests.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Mice Tested")
plt.show()


# In[27]:


# Creating the same bar plot using pandas + Pyplot.

# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.
tests = no_dupes_df["Drug Regimen"].value_counts()
plt.bar(tests.index.values, tests.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Mice Tested")
plt.show()


# In[28]:


# Creating a pie plot using pandas dataframe.

# Generate a pie plot showing the distribution of female versus male mice using Pandas
gender = no_dupes_df.Sex.value_counts()
gender.plot(kind="pie", autopct="%1.1f%%")


# In[29]:


# Creating the same pie plot using pandas + Pyplot.

# Generate a pie plot showing the distribution of female versus male mice using pyplot
gender = no_dupes_df.Sex.value_counts()
plt.pie(gender.values, labels=gender.index.values, autopct="%1.1f%%")
plt.ylabel("Sex")
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[35]:


# Starting to analyse data for outliers in select treatments.

# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
last_tumor=no_dupes_df.groupby(["Mouse ID"])["Timepoint"].max()
last_tumor=last_tumor.reset_index()
# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
tumor_merged=last_tumor.merge(no_dupes_df, on=["Mouse ID", "Timepoint"],how="left")


# In[40]:


# Detecting potential outliers in the treatment, along with organizing data for plots later.

# Put treatments into a list for for loop (and later for plot labels)
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_list = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for drug in treatment_list:
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    last_tumor_vol = tumor_merged.loc[tumor_merged["Drug Regimen"]==drug,"Tumor Volume (mm3)"]
    
    # add subset 
    tumor_vol_list.append(last_tumor_vol)
    
    # Determine outliers using upper and lower bounds
    quartiles = last_tumor_vol.quantile([.25,.5,.75])
    firstq = quartiles[0.25]
    thirdq = quartiles[0.75]
    iqr = thirdq - firstq
    lower_bound = firstq - (1.5 * iqr)
    upper_bound = thirdq + (1.5 * iqr)
    
    outliers = last_tumor_vol.loc[(last_tumor_vol < lower_bound) | (last_tumor_vol > upper_bound)]
    print(f"{drug}'s potential outliers: {outliers}")


# In[42]:


# Creating the box plot with visible outliers based on data from prior code.

# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
orange_out = dict(markerfacecolor='red', markersize=10)
plt.boxplot(tumor_vol_list, labels = treatment_list, flierprops=orange_out)
plt.ylabel("Last Tumor Volume (mm3)")
plt.show()


# ## Line and Scatter Plots

# In[50]:


# Collecting data for mouse L509 for plots like this line plot.

# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
capomulin_stats = no_dupes_df[no_dupes_df["Drug Regimen"] == "Capomulin"]
mouse_data = capomulin_stats[capomulin_stats["Mouse ID"] == "l509"]
plt.plot(mouse_data["Timepoint"], mouse_data["Tumor Volume (mm3)"])
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.title("Capomulin Treatment of Mouse l509")


# In[54]:


# Collecting data for Capomulin treated mice for plots like this scatter plot.

# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
capomulin_stats = no_dupes_df[no_dupes_df["Drug Regimen"] == "Capomulin"]
capomulin_average = capomulin_stats.groupby(["Mouse ID"]).mean()
plt.scatter(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.title("Tumor Volume vs Weight of Capomulin Treated Mice")


# ## Correlation and Regression

# In[59]:


# Adding a linear regression model into recently created scatter plot.

# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
corre = st.pearsonr(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
print(f"The correlation between mouse weight and the average tumor volume is {round(corre[0],2)}")

model = st.linregress(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
slope = model [0]
y_intercept = model [1]
y_values = capomulin_average["Weight (g)"] * slope + y_intercept
plt.scatter(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.title("Tumor Volume vs Weight of Capomulin Treated Mice")
plt.plot(capomulin_average["Weight (g)"], y_values, color="red")
