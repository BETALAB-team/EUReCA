import os
import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.style.use('ggplot')
from sklearn.metrics import mean_squared_error
# ggplot, bmh, seaborn, classic, default



tests_dict = pd.read_excel(os.path.join(".","validation_Test_VDI.xlsx"), sheet_name=None)
TESTS = [1,3,
         2,4,
         5,
         6,
         7,
         12
    ]
RMSE = pd.DataFrame(index = [f"Test {t}" for t in TESTS], columns = pd.MultiIndex.from_product([["T","Load"],["Day 1","Day 10","Day 60"]]))
for test_id in TESTS:
    columns = pd.MultiIndex.from_product([["Day 1","Day 10","Day 60"],["T 1C","T 2C","Top 1C","Top 2C", "T 6007", "Top 6007", "T 6020", "Load 1C","Load 2C", "Load 6007","Load 6020"]])
    data = pd.DataFrame(0., index = range(24), columns = columns)


    results = pd.read_csv(os.path.join("VDI_tests_results",f"Tests_{test_id}.csv"), header = [0,1], index_col = [0])

    data["Day 1","T 1C"] = results["1C"]["Ta [°C]"].iloc[0:24].values
    data["Day 1", "T 2C"] = results["2C"]["Ta [°C]"].iloc[0:24].values
    data["Day 10","T 1C"] = results["1C"]["Ta [°C]"].iloc[24*9:24*10].values
    data["Day 10", "T 2C"] = results["2C"]["Ta [°C]"].iloc[24*9:24*10].values
    data["Day 60","T 1C"] = results["1C"]["Ta [°C]"].iloc[24*59:24*60].values
    data["Day 60", "T 2C"] = results["2C"]["Ta [°C]"].iloc[24*59:24*60].values
    data["Day 1", "Top 1C"] = results["1C"]["Top [°C]"].iloc[0:24].values
    data["Day 1", "Top 2C"] = results["2C"]["Top [°C]"].iloc[0:24].values
    data["Day 10", "Top 1C"] = results["1C"]["Top [°C]"].iloc[24 * 9:24 * 10].values
    data["Day 10", "Top 2C"] = results["2C"]["Top [°C]"].iloc[24 * 9:24 * 10].values
    data["Day 60", "Top 1C"] = results["1C"]["Top [°C]"].iloc[24 * 59:24 * 60].values
    data["Day 60", "Top 2C"] = results["2C"]["Top [°C]"].iloc[24 * 59:24 * 60].values
    data["Day 1","Load 1C"] = results["1C"]["Sens Load [W]"].iloc[0:24].values
    data["Day 1", "Load 2C"] = results["2C"]["Sens Load [W]"].iloc[0:24].values
    data["Day 10","Load 1C"] = results["1C"]["Sens Load [W]"].iloc[24*9:24*10].values
    data["Day 10", "Load 2C"] = results["2C"]["Sens Load [W]"].iloc[24*9:24*10].values
    data["Day 60","Load 1C"] = results["1C"]["Sens Load [W]"].iloc[24*59:24*60].values
    data["Day 60", "Load 2C"] = results["2C"]["Sens Load [W]"].iloc[24*59:24*60].values

    data["Day 1", "T 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 6"].values
    data["Day 1", "T 6020"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 1"].values
    data["Day 10", "T 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 7"].values
    data["Day 10", "T 6020"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 2"].values
    data["Day 60", "T 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 8"].values
    data["Day 60", "T 6020"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 3"].values
    data["Day 1", "Top 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 9"].values
    data["Day 10", "Top 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 10"].values
    data["Day 60", "Top 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 11"].values
    data["Day 1", "Load 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 19"].values
    data["Day 1", "Load 6020"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 13"].values
    data["Day 10", "Load 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 20"].values
    data["Day 10", "Load 6020"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 14"].values
    data["Day 60", "Load 6007"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 21"].values
    data["Day 60", "Load 6020"] = tests_dict[f"Test{test_id}"].iloc[3:27]["Unnamed: 15"].values

    RMSE.loc[f"Test {test_id}"]["T","Day 1"]  = mean_squared_error(data["Day 1", "T 6007"],data["Day 1", "T 2C"], squared = False)
    RMSE.loc[f"Test {test_id}"]["T","Day 10"]  = mean_squared_error(data["Day 10", "T 6007"],data["Day 10", "T 2C"], squared = False)
    RMSE.loc[f"Test {test_id}"]["T","Day 60"]  = mean_squared_error(data["Day 60", "T 6007"],data["Day 60", "T 2C"], squared = False)
    RMSE.loc[f"Test {test_id}"]["Load","Day 1"] = mean_squared_error(data["Day 1", "Load 6007"], data["Day 1", "Load 2C"],
                                                                  squared=False)
    RMSE.loc[f"Test {test_id}"]["Load","Day 10"] = mean_squared_error(data["Day 10", "Load 6007"], data["Day 10", "Load 2C"],
                                                                   squared=False)
    RMSE.loc[f"Test {test_id}"]["Load","Day 60"] = mean_squared_error(data["Day 60", "Load 6007"], data["Day 60", "Load 2C"],
                                                                   squared=False)

    fig, [[ax11,ax12],[ax21,ax22],[ax31,ax32]] = plt.subplots(nrows=3, ncols=2, figsize = (12,12))
    style = [':','-','-','-', "-.", "-."]
    print_t = ["T 1C","T 2C", "T 6007","T 6020", "Top 2C", "Top 6007"]
    print_load = ["Load 1C","Load 2C", "Load 6007","Load 6020"]
    data["Day 1"][print_t].plot(ax = ax11, style = style)
    data["Day 10"][print_t].plot(ax = ax21, style = style)
    data["Day 60"][print_t].plot(ax = ax31, style = style)
    data["Day 1"][print_load].plot(ax = ax12, style = style)
    data["Day 10"][print_load].plot(ax = ax22, style = style)
    data["Day 60"][print_load].plot(ax = ax32, style = style)

    fig.suptitle(f"VDI 6007 BESTEST: TEST {test_id}")


    for ax, day in zip([ax11,ax21,ax31],[1,10,60]):
        average = data[f"Day {day}"][print_t].mean().mean()
        ax.set_title(f"Day {day}")
        ax.set_ylabel(f"Temperature [°C]")
        ax.grid(alpha = 0.5)
        ax.set_xticks(ax.get_xticks().tolist())
        ax.set_xticklabels([f"{int(h)}:00" for h in ax.get_xticks().tolist()])
        ax.set_xlim([-.5,24.5])
        ax.set_ylim([average-6.5,average+6.5])


    for ax, day in zip([ax12,ax22,ax32],[1,10,60]):
        ax.set_title(f"Day {day}")
        ax.set_ylabel(f"Load [W]")
        ax.grid(alpha = 0.5)
        ax.set_xticks(ax.get_xticks().tolist())
        ax.set_xticklabels([f"{int(h)}:00" for h in ax.get_xticks().tolist()])
        ax.set_xlim([-.5,24.5])

    plt.tight_layout()

    fig.savefig(os.path.join("VDI_tests_results",f"Test_{test_id}.png"))
    plt.close()

RMSE["T"].plot(kind = "bar")
RMSE["Load"].plot(kind = "bar")