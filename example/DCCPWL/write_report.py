
def solving_time_synthesis(report_file):
    # read the report file
    time={
        "MA": 0,
        "CPWL-MA": 0,
        "QCPWL-MA": 0,
        "Subgrad": 0,
        "MINLP": 0,
    }
    with open(report_file, "r") as logfile:
        algo_name = None
        for line in logfile:
            if line.startswith("Testing MA"):
                algo_name = "MA"
            elif line.startswith("Testing CPWL-MA"):
                algo_name = "CPWL-MA"
            elif line.startswith("Testing QCPWL-MA"):
                algo_name = "QCPWL-MA"
            elif line.startswith("Testing QCPWL-Subgrad-MINLP"):
                algo_name = "MINLP"
            elif line.startswith("Testing QCPWL-Subgrad"):
                algo_name = "Subgrad"
            elif line.startswith("Testing QCPWL-MINLP"):
                algo_name = "MINLP"
            if line.startswith("Total seconds in IPOPT"):
                time[algo_name] += float(line.split()[-1])
            elif line.startswith("Deterministic time = "):
                line_split = line.split()
                ticks = float(line_split[3])
                ticks_per_sec = float(line_split[5].replace("(", ""))
                time[algo_name] += ticks/ticks_per_sec
            elif line.startswith("Solving Time (sec)"):
                time[algo_name] += float(line.split()[-1])
    print(time)


if __name__ == "__main__":
    solving_time_synthesis(r"report/unconstraint_wo_for paper 20230524.txt")
    solving_time_synthesis(r"report/HPC_for paper 20230524.txt")
