#include "Solver.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>

#include <cmath>

#include "MpSolver.h"


using namespace std;


namespace szx {

#pragma region Solver::Cli
int Solver::Cli::run(int argc, char * argv[]) {
    Log(LogSwitch::Szx::Cli) << "parse command line arguments." << endl;
    Set<String> switchSet;
    Map<String, char*> optionMap({ // use string as key to compare string contents instead of pointers.
        { InstancePathOption(), nullptr },
        { SolutionPathOption(), nullptr },
        { RandSeedOption(), nullptr },
        { TimeoutOption(), nullptr },
        { MaxIterOption(), nullptr },
        { JobNumOption(), nullptr },
        { RunIdOption(), nullptr },
        { EnvironmentPathOption(), nullptr },
        { ConfigPathOption(), nullptr },
        { LogPathOption(), nullptr }
    });

    for (int i = 1; i < argc; ++i) { // skip executable name.
        auto mapIter = optionMap.find(argv[i]);
        if (mapIter != optionMap.end()) { // option argument.
            mapIter->second = argv[++i];
        } else { // switch argument.
            switchSet.insert(argv[i]);
        }
    }

    Log(LogSwitch::Szx::Cli) << "execute commands." << endl;
    if (switchSet.find(HelpSwitch()) != switchSet.end()) {
        cout << HelpInfo() << endl;
    }

    if (switchSet.find(AuthorNameSwitch()) != switchSet.end()) {
        cout << AuthorName() << endl;
    }

    Solver::Environment env;
    env.load(optionMap);
    if (env.instPath.empty() || env.slnPath.empty()) { return -1; }

    Solver::Configuration cfg;
    cfg.load(env.cfgPath);

    Log(LogSwitch::Szx::Input) << "load instance " << env.instPath << " (seed=" << env.randSeed << ")." << endl;
    Problem::Input input;
    if (!input.load(env.instPath)) { return -1; }

    Solver solver(input, env, cfg);
    solver.solve();

    pb::Submission submission;
    submission.set_thread(to_string(env.jobNum));
    submission.set_instance(env.friendlyInstName());
    submission.set_duration(to_string(solver.timer.elapsedSeconds()) + "s");

    solver.output.save(env.slnPath, submission);
    #if SZX_DEBUG
    solver.output.save(env.solutionPathWithTime(), submission);
    solver.record();
    #endif // SZX_DEBUG

    return 0;
}
#pragma endregion Solver::Cli

#pragma region Solver::Environment
void Solver::Environment::load(const Map<String, char*> &optionMap) {
    char *str;

    str = optionMap.at(Cli::EnvironmentPathOption());
    if (str != nullptr) { loadWithoutCalibrate(str); }

    str = optionMap.at(Cli::InstancePathOption());
    if (str != nullptr) { instPath = str; }

    str = optionMap.at(Cli::SolutionPathOption());
    if (str != nullptr) { slnPath = str; }

    str = optionMap.at(Cli::RandSeedOption());
    if (str != nullptr) { randSeed = atoi(str); }

    str = optionMap.at(Cli::TimeoutOption());
    if (str != nullptr) { msTimeout = static_cast<Duration>(atof(str) * Timer::MillisecondsPerSecond); }

    str = optionMap.at(Cli::MaxIterOption());
    if (str != nullptr) { maxIter = atoi(str); }

    str = optionMap.at(Cli::JobNumOption());
    if (str != nullptr) { jobNum = atoi(str); }

    str = optionMap.at(Cli::RunIdOption());
    if (str != nullptr) { rid = str; }

    str = optionMap.at(Cli::ConfigPathOption());
    if (str != nullptr) { cfgPath = str; }

    str = optionMap.at(Cli::LogPathOption());
    if (str != nullptr) { logPath = str; }

    calibrate();
}

void Solver::Environment::load(const String &filePath) {
    loadWithoutCalibrate(filePath);
    calibrate();
}

void Solver::Environment::loadWithoutCalibrate(const String &filePath) {
    // EXTEND[szx][8]: load environment from file.
    // EXTEND[szx][8]: check file existence first.
}

void Solver::Environment::save(const String &filePath) const {
    // EXTEND[szx][8]: save environment to file.
}
void Solver::Environment::calibrate() {
    // adjust thread number.
    int threadNum = thread::hardware_concurrency();
    if ((jobNum <= 0) || (jobNum > threadNum)) { jobNum = threadNum; }

    // adjust timeout.
    msTimeout -= Environment::SaveSolutionTimeInMillisecond;
}
#pragma endregion Solver::Environment

#pragma region Solver::Configuration
void Solver::Configuration::load(const String &filePath) {
    // EXTEND[szx][5]: load configuration from file.
    // EXTEND[szx][8]: check file existence first.
}

void Solver::Configuration::save(const String &filePath) const {
    // EXTEND[szx][5]: save configuration to file.
}
#pragma endregion Solver::Configuration

#pragma region Solver
bool Solver::solve() {
    init();

    int workerNum = (max)(1, env.jobNum / cfg.threadNumPerWorker);
    cfg.threadNumPerWorker = env.jobNum / workerNum;
    List<Solution> solutions(workerNum, Solution(this));
    List<bool> success(workerNum);

    Log(LogSwitch::Szx::Framework) << "launch " << workerNum << " workers." << endl;
    List<thread> threadList;
    threadList.reserve(workerNum);
    for (int i = 0; i < workerNum; ++i) {
        // TODO[szx][2]: as *this is captured by ref, the solver should support concurrency itself, i.e., data members should be read-only or independent for each worker.
        // OPTIMIZE[szx][3]: add a list to specify a series of algorithm to be used by each threads in sequence.
        threadList.emplace_back([&, i]() { success[i] = optimize(solutions[i], i); });
    }
    for (int i = 0; i < workerNum; ++i) { threadList.at(i).join(); }

    Log(LogSwitch::Szx::Framework) << "collect best result among all workers." << endl;
    int bestIndex = -1;
    Length bestValue = Problem::MaxColorNum;
    for (int i = 0; i < workerNum; ++i) {
        if (!success[i]) { continue; }
        Log(LogSwitch::Szx::Framework) << "worker " << i << " got " << solutions[i].colorNum << endl;
        if (solutions[i].colorNum >= bestValue) { continue; }
        bestIndex = i;
        bestValue = solutions[i].colorNum;
    }

    env.rid = to_string(bestIndex);
    if (bestIndex < 0) { return false; }
    output = solutions[bestIndex];
    return true;
}

void Solver::record() const {
    #if SZX_DEBUG
    int generation = 0;

    ostringstream log;

    System::MemoryUsage mu = System::peakMemoryUsage();

    Length obj = output.colorNum;
    Length checkerObj = -1;
    bool feasible = check(checkerObj);

    // record basic information.
    log << env.friendlyLocalTime() << ","
        << env.rid << ","
        << env.instPath << ","
        << feasible << "," << (obj - checkerObj) << ","
        << obj << ","
        << timer.elapsedSeconds() << ","
        << mu.physicalMemory << "," << mu.virtualMemory << ","
        << env.randSeed << ","
        << cfg.toBriefStr() << ","
        << generation << "," << iteration << ",";

    // record solution vector.
    for (auto n = output.nodecolors().begin(); n != output.nodecolors().end(); ++n) {
        log << *n << " ";
    }
    log << endl;

    // append all text atomically.
    static mutex logFileMutex;
    lock_guard<mutex> logFileGuard(logFileMutex);

    ofstream logFile(env.logPath, ios::app);
    logFile.seekp(0, ios::end);
    if (logFile.tellp() <= 0) {
        logFile << "Time,ID,Instance,Feasible,ObjMatch,Color,Duration,PhysMem,VirtMem,RandSeed,Config,Generation,Iteration,Solution" << endl;
    }
    logFile << log.str();
    logFile.close();
    #endif // SZX_DEBUG
}

bool Solver::check(Length &checkerObj) const {
    #if SZX_DEBUG
    enum CheckerFlag {
        IoError = 0x0,
        FormatError = 0x1,
        ColorConflictError = 0x2
    };

    checkerObj = System::exec("Checker.exe " + env.instPath + " " + env.solutionPathWithTime());
    if (checkerObj > 0) { return true; }
    checkerObj = ~checkerObj;
    if (checkerObj == CheckerFlag::IoError) { Log(LogSwitch::Checker) << "IoError." << endl; }
    if (checkerObj & CheckerFlag::FormatError) { Log(LogSwitch::Checker) << "FormatError." << endl; }
    if (checkerObj & CheckerFlag::ColorConflictError) { Log(LogSwitch::Checker) << "ColorConflictError." << endl; }
    return false;
    #else
    checkerObj = 0;
    return true;
    #endif // SZX_DEBUG
}

void Solver::init() {
    ID nodeNum = input.graph().nodenum();

    aux.adjList.init(nodeNum);
    for (auto e = input.graph().edges().begin(); e != input.graph().edges().end(); ++e) {
        // assume there is no duplicated edge.
        aux.adjList.at(e->src()).push_back(e->dst());
        aux.adjList.at(e->src()).push_back(e->dst());
    }
}

bool Solver::optimize(Solution &sln, ID workerId) {
    Log(LogSwitch::Szx::Framework) << "worker " << workerId << " starts." << endl;

    // reset solution state.
    bool status = false;

    //status = optimizeBoolDecisionModel(sln);
    //status = optimizeRelaxedBoolDecisionModel(sln);
    status = optimizeIntegerDecisionModel(sln);

    sln.colorNum = input.colornum(); // record obj.

    Log(LogSwitch::Szx::Framework) << "worker " << workerId << " ends." << endl;
    return status;
}

bool Solver::optimizeBoolDecisionModel(Solution &sln) {
    ID nodeNum = input.graph().nodenum();

    auto &nodeColors(*sln.mutable_nodecolors());
    nodeColors.Resize(nodeNum, Problem::InvalidId);

    MpSolver mp;

    // add decision variables.
    Arr2D<MpSolver::DecisionVar> isColor(nodeNum, input.colornum());
    for (auto n = 0; n < nodeNum; ++n) {
        for (int c = 0; c < input.colornum(); ++c) {
            isColor.at(n, c) = mp.addVar(MpSolver::VariableType::Bool, 0, 1, 0);
        }
    }

    // add constraints.
    // single color.
    for (ID n = 0; n < nodeNum; ++n) {
        MpSolver::LinearExpr sum;
        for (auto c = 0; c < input.colornum(); ++c) {
            sum += isColor.at(n, c);
        }
        //mp.addConstraint(sum >= 1);
        mp.addConstraint(sum == 1);
    }

    // conflict avoidance.
    for (ID n = 0; n < nodeNum; ++n) {
        for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
            for (auto c = 0; c < input.colornum(); ++c) {
                mp.addConstraint(isColor.at(n, c) + isColor.at(*m, c) <= 1);
            }
        }
    }

    // solve model.
    mp.setOutput(true);
    //mp.setMaxThread(1);
    mp.setTimeLimitInSecond(1800);
    //mp.setMipFocus(MpSolver::MipFocusMode::ImproveFeasibleSolution);

    // record decision.
    if (mp.optimize()) {
        for (ID n = 0; n < nodeNum; ++n) {
            for (ID c = 0; c < input.colornum(); ++c) {
                if (mp.isTrue(isColor.at(n, c))) { nodeColors[n] = c; break; }
            }
        }
        return true;
    }

    return false;
}

bool Solver::optimizeRelaxedBoolDecisionModel(Solution &sln) {
    ID nodeNum = input.graph().nodenum();

    auto &nodeColors(*sln.mutable_nodecolors());
    nodeColors.Resize(nodeNum, Problem::InvalidId);

    MpSolver mp;

    // add decision variables.
    Arr2D<MpSolver::DecisionVar> isColor(nodeNum, input.colornum());
    for (auto n = 0; n < nodeNum; ++n) {
        for (int c = 0; c < input.colornum(); ++c) {
            isColor.at(n, c) = mp.addVar(MpSolver::VariableType::Bool, 0, 1, 0);
        }
    }
    Arr2D<MpSolver::DecisionVar> hasConflict(nodeNum, nodeNum);
    for (auto n = 0; n < nodeNum; ++n) {
        for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
            hasConflict.at(n, *m) = mp.addVar(MpSolver::VariableType::Real, 0, 1, 1);
        }
    }

    // add constraints.
    // single color.
    for (ID n = 0; n < nodeNum; ++n) {
        MpSolver::LinearExpr sum;
        for (auto c = 0; c < input.colornum(); ++c) {
            sum += isColor.at(n, c);
        }
        //mp.addConstraint(sum >= 1);
        mp.addConstraint(sum == 1);
    }

    // conflict avoidance.
    for (ID n = 0; n < nodeNum; ++n) {
        for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
            for (auto c = 0; c < input.colornum(); ++c) {
                mp.addConstraint(isColor.at(n, c) + isColor.at(*m, c) <= 1 + hasConflict.at(n, *m));
            }
        }
    }

    // set objective.
    mp.setOptimaOrientation(MpSolver::OptimaOrientation::Minimize);

    // solve model.
    mp.setOutput(true);
    //mp.setMaxThread(1);
    mp.setTimeLimitInSecond(1800);
    //mp.setMipFocus(MpSolver::MipFocusMode::ImproveFeasibleSolution);

    // record decision.
    if (mp.optimize()) {
        for (ID n = 0; n < nodeNum; ++n) {
            for (ID c = 0; c < input.colornum(); ++c) {
                if (mp.isTrue(isColor.at(n, c))) { nodeColors[n] = c; break; }
            }
        }
        return true;
    }

    return false;
}

bool Solver::optimizeIntegerDecisionModel(Solution &sln) {
    ID nodeNum = input.graph().nodenum();

    auto &nodeColors(*sln.mutable_nodecolors());
    nodeColors.Resize(nodeNum, Problem::InvalidId);

    MpSolver mp;

    // add decision variables.
    Arr<MpSolver::DecisionVar> colors(nodeNum);
    for (auto x = colors.begin(); x != colors.end(); ++x) {
        *x = mp.addVar(MpSolver::VariableType::Integer, 0, input.colornum() - 1, 0);
    }
    Arr2D<MpSolver::DecisionVar> isSmallerColor(nodeNum, nodeNum);
    for (auto n = 0; n < nodeNum; ++n) {
        for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
            isSmallerColor.at(n, *m) = mp.addVar(MpSolver::VariableType::Bool, 0, 1, 0);
        }
    }

    // add constraints.
    // conflict avoidance.
    for (ID n = 0; n < nodeNum; ++n) {
        for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
            MpSolver::LinearExpr diff = colors[n] - colors[*m];
            mp.addConstraint(1 - input.colornum() * (1 - isSmallerColor.at(n, *m)) <= diff);
            mp.addConstraint(diff <= input.colornum() * isSmallerColor.at(n, *m) - 1);
        }
    }

    // solve model.
    mp.setOutput(true);
    //mp.setMaxThread(1);
    mp.setTimeLimitInSecond(1800);
    //mp.setMipFocus(MpSolver::MipFocusMode::ImproveFeasibleSolution);

    // record decision.
    if (mp.optimize()) {
        for (ID n = 0; n < nodeNum; ++n) {
            nodeColors[n] = lround(mp.getValue(colors.at(n)));
        }
        return true;
    }

    return false;
}

void Solver::detectClique() {
    aux.fixedColors.resize(input.graph().nodenum(), Problem::InvalidId);

    auto e = input.graph().edges().begin();
    tsm::solveWeightedMaxClique(aux.clique, [&](ID &src, ID &dst) {
        if (e == input.graph().edges().end()) {
            e = input.graph().edges().begin();
            return false;
        }
        src = e->src();
        dst = e->dst();
        ++e;
        return true;
    }, input.graph().nodenum());

    Log(LogSwitch::Szx::Preprocess) << "clique[" << aux.clique.weight << "]=";
    ID c = 0;
    for (auto n = aux.clique.nodes.begin(); n != aux.clique.nodes.end(); ++n) {
        aux.fixedColors[*n] = c++;
        Log(LogSwitch::Szx::Preprocess) << " " << *n;
    }
    Log(LogSwitch::Szx::Preprocess) << endl;
}
#pragma endregion Solver

}
