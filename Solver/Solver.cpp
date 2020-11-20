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
    aux.adjMat.init(nodeNum, nodeNum);
    aux.adjMat.reset(Arr2D<bool>::ResetOption::AllBits0);
    for (auto e = input.graph().edges().begin(); e != input.graph().edges().end(); ++e) {
        // assume there is no duplicated edge.
        aux.adjList.at(e->src()).push_back(e->dst());
        aux.adjList.at(e->dst()).push_back(e->src());
        aux.adjMat.at(e->src(), e->dst()) = true;
        aux.adjMat.at(e->dst(), e->src()) = true;
    }

    detectClique();
}

bool Solver::optimize(Solution &sln, ID workerId) {
    Log(LogSwitch::Szx::Framework) << "worker " << workerId << " starts." << endl;

    // reset solution state.
    bool status = false;

    //status = optimizeBoolDecisionModel(sln);
    //status = optimizeRelaxedBoolDecisionModel(sln);
    //status = optimizeIntegerDecisionModel(sln);
    //status = optimizeCoveringRelaxedBoolDecisionModel(sln);
    //status = optimizeCoveringBoolDecisionModel(sln);
    //status = optimizeLocalSearch(sln);
    status = optimizeTabuSearchPQ(sln);
    //status = optimizeTabuSearch(sln);
    //status = optimizeAuctionSearch(sln);
    //status = optimizeCliqueReduction(sln);

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
    //detectIndependentSetsInSubGraphs();

    ID nodeNum = input.graph().nodenum();

    auto &nodeColors(*sln.mutable_nodecolors());
    nodeColors.Resize(nodeNum, Problem::InvalidId);

    MpSolver mp;

    // add decision variables.
    Arr2D<MpSolver::DecisionVar> isColor(nodeNum, input.colornum());
    for (auto n = 0; n < nodeNum; ++n) {
        bool isFixed = aux.fixedColors[n] > Problem::InvalidId;
        for (int c = 0; c < input.colornum(); ++c) {
            double lb = isFixed ? (aux.fixedColors[n] == c) : 0;
            double ub = isFixed ? (aux.fixedColors[n] == c) : 1;
            isColor.at(n, c) = mp.addVar(MpSolver::VariableType::Bool, lb, ub, 0);
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
    //for (ID n = 0; n < nodeNum; ++n) {
    //    ID maxSameColorNum = static_cast<ID>(aux.independentSetInSubGraphs[n].nodes.size());
    //    for (auto c = 0; c < input.colornum(); ++c) {
    //        MpSolver::LinearExpr degree;
    //        for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
    //            degree += isColor.at(*m, c);
    //        }
    //        mp.addConstraint(maxSameColorNum * isColor.at(n, c) + degree <= maxSameColorNum + maxSameColorNum * hasConflict.at(n, aux.adjList[n].front()));
    //    }
    //}

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

bool Solver::optimizeCoveringBoolDecisionModel(Solution &sln) {
    detectIndependentSetsInInvSubGraphs();

    ID nodeNum = input.graph().nodenum();

    Arr<double> coverWeights(nodeNum, 1); // for node coverage.
    Arr<tsm::Weight> iSetWeights(nodeNum, 1); // for independent set computing.

    while (!optimizeCoveringBoolDecisionModel(sln, coverWeights, iSetWeights, true)) {
        //for (ID n = 0; n < nodeNum; ++n) { // increase weights for uncovered nodes.
        //    //if (sln.nodecolors(n) >= 0) { continue; }
        //    //iSetWeights[n] *= 8;
        //    //coverWeights[n] += 0.125;
        //    if (sln.nodecolors(n) < 0) {
        //        iSetWeights[n] = 16;
        //        coverWeights[n] = 1.0625;
        //        Log(LogSwitch::Szx::Model) << " " << n;
        //    } else {
        //        iSetWeights[n] = 1;
        //        coverWeights[n] = 1;
        //    }
        //}
        //Log(LogSwitch::Szx::Model) << endl;

        tsm::Clique iSet;
        tsm::solveWeightedIndependentSet(iSet, aux.adjMat, iSetWeights, cfg.msCliqueDetectionTimeout);

        tsm::Weight error = static_cast<tsm::Weight>(iSet.nodes.size());
        tsm::Weight reducedCost = iSet.weight - error;
        Log(LogSwitch::Szx::Model) << "ReducedCost=" << reducedCost << endl;

        sort(iSet.nodes.begin(), iSet.nodes.end());
        if (aux.independentSets.set(iSet.nodes, iSet)) {
            Log(LogSwitch::Szx::Model) << "IndependentSet[" << iSet.nodes.size() << "]=";
            for (auto m = iSet.nodes.begin(); m != iSet.nodes.end(); ++m) {
                Log(LogSwitch::Szx::Model) << " " << *m;
            }
            Log(LogSwitch::Szx::Model) << endl;
        } else if (reducedCost <= Configuration::TsmPrecision + (error / 2)) {
            break;
        }
    }

    optimizeCoveringRelaxedBoolDecisionModel(sln, coverWeights, iSetWeights);

    return true;
}

bool Solver::optimizeCoveringBoolDecisionModel(Solution &sln, Arr<double> &coverWeights, Arr<tsm::Weight> &iSetWeights, bool linearRelax) {
    detectIndependentSetsInInvSubGraphs();

    ID nodeNum = input.graph().nodenum();

    auto &nodeColors(*sln.mutable_nodecolors());
    nodeColors.Clear();
    nodeColors.Resize(nodeNum, Problem::InvalidId);

    Arr<List<ID>> nodeInSets(nodeNum); // nodeInSets[n] is a list of independent sets containing node n.
    ID iSetId = 0;
    aux.independentSets.forEach([&](const tsm::Clique &iSet) {
        for (auto n = iSet.nodes.begin(); n != iSet.nodes.end(); ++n) {
            nodeInSets[*n].push_back(iSetId);
        }
        ++iSetId;
        return false;
    });

    MpSolver mp;

    // add decision variables.
    MpSolver::VariableType varType = linearRelax ? MpSolver::VariableType::Real : MpSolver::VariableType::Bool;
    Arr<MpSolver::DecisionVar> isPicked(aux.independentSets.dataNum());
    for (ID s = 0; s < aux.independentSets.dataNum(); ++s) {
        isPicked[s] = mp.addVar(varType, 0, 1, 1);
    }

    // add constraints.
    Arr<MpSolver::Constraint> coverages(nodeNum);
    MpSolver::Constraint setNum;

    // node coverage.
    for (ID n = 0; n < nodeNum; ++n) {
        MpSolver::LinearExpr coveringSetNum;
        for (auto s = nodeInSets[n].begin(); s != nodeInSets[n].end(); ++s) {
            coveringSetNum += isPicked[*s];
        }
        coverages[n] = mp.addConstraint(coveringSetNum >= 1);
    }

    //// set number.
    //MpSolver::LinearExpr pickedSetNum;
    //for (ID s = 0; s < aux.independentSets.dataNum(); ++s) {
    //    pickedSetNum += isPicked[s];
    //}
    ////mp.addConstraint(pickedSetNum == input.colornum());
    //setNum = mp.addConstraint(pickedSetNum <= input.colornum());

    // set objective.
    mp.setOptimaOrientation(MpSolver::OptimaOrientation::Minimize);

    // solve model.
    mp.setOutput(true);
    //mp.setMaxThread(1);
    mp.setTimeLimitInSecond(1800);
    //mp.setMipFocus(MpSolver::MipFocusMode::ImproveFeasibleSolution);

    // record decision.
    if (mp.optimize()) {
        if (linearRelax) { // record shadow price.
            for (ID n = 0; n < nodeNum; ++n) {
                iSetWeights[n] = lround(Configuration::TsmPrecision * mp.getDualValue(coverages[n])) + 1;
            }
        } else { // collect integer solution.
            ID coloredNodeNum = 0;
            ID color = 0;
            ID s = 0;
            aux.independentSets.forEach([&](const tsm::Clique &iSet) {
                if (!mp.isTrue(isPicked[s++])) { return false; }
                for (auto n = iSet.nodes.begin(); n != iSet.nodes.end(); ++n) {
                    if (nodeColors[*n] >= 0) { continue; }
                    nodeColors[*n] = color;
                    ++coloredNodeNum;
                }
                return (++color >= input.colornum()); // stop assigning colors if the colors are used up.
            });
            if (coloredNodeNum >= nodeNum) { return true; }
            Log(LogSwitch::Szx::Model) << "Uncovered[" << (nodeNum - coloredNodeNum) << "]=";
            for (ID n = 0; n < nodeNum; ++n) {
                if (sln.nodecolors(n) < 0) { Log(LogSwitch::Szx::Model) << " " << n; }
            }
            Log(LogSwitch::Szx::Model) << endl;
        }

    }

    return false;
}

bool Solver::optimizeCoveringRelaxedBoolDecisionModel(Solution &sln) {
    detectIndependentSetsInInvSubGraphs();

    ID nodeNum = input.graph().nodenum();

    Arr<double> coverWeights(nodeNum, 1); // for node coverage.
    Arr<tsm::Weight> iSetWeights(nodeNum, 1); // for independent set computing.

    while (!optimizeCoveringRelaxedBoolDecisionModel(sln, coverWeights, iSetWeights)) {
        for (ID n = 0; n < nodeNum; ++n) { // increase weights for uncovered nodes.
            //iSetWeights[n] *= 8;
            //coverWeights[n] += 0.125;
            if (sln.nodecolors(n) < 0) {
                iSetWeights[n] = 16;
                coverWeights[n] = 1.0625;
                Log(LogSwitch::Szx::Model) << " " << n;
            } else {
                iSetWeights[n] = 1;
                coverWeights[n] = 1;
            }
        }
        Log(LogSwitch::Szx::Model) << endl;

        tsm::Clique iSet;
        tsm::solveWeightedIndependentSet(iSet, aux.adjMat, iSetWeights, cfg.msCliqueDetectionTimeout);

        sort(iSet.nodes.begin(), iSet.nodes.end());
        if (aux.independentSets.set(iSet.nodes, iSet)) {
            Log(LogSwitch::Szx::Model) << "IndependentSet[" << iSet.nodes.size() << "]=";
            for (auto m = iSet.nodes.begin(); m != iSet.nodes.end(); ++m) {
                Log(LogSwitch::Szx::Model) << " " << *m;
            }
            Log(LogSwitch::Szx::Model) << endl;
        }
    }

    //while (!optimizeCoveringRelaxedBoolDecisionModel(sln, coverWeights, iSetWeights)) {
    //    for (ID n = 0; n < nodeNum; ++n) { // increase weights for uncovered nodes.
    //        if (sln.nodecolors(n) >= 0) { continue; }
    //        iSetWeights[n] *= 8;
    //        coverWeights[n] += 0.125;
    //        Log(LogSwitch::Szx::Model) << " " << n;
    //    }
    //    Log(LogSwitch::Szx::Model) << endl;

    //    tsm::Clique iSet;
    //    for (ID n = 0; n < nodeNum; ++n) {
    //        Subgraph &subInvGraph(aux.subInvGraphs[n]);
    //        Arr<tsm::Weight> weights(static_cast<ID>(subInvGraph.idMap.size()));
    //        for (ID i = 0; i < weights.size(); ++i) { weights[i] = iSetWeights[subInvGraph.idMap[i]]; }

    //        tsm::solveWeightedMaxClique(iSet, subInvGraph.adjMat, weights, cfg.msCliqueDetectionTimeout);

    //        subInvGraph.mapBack(iSet.nodes);
    //        sort(iSet.nodes.begin(), iSet.nodes.end());
    //        if (aux.independentSets.set(iSet.nodes, iSet)) {
    //            Log(LogSwitch::Szx::Model) << n << ".IndependentSet[" << iSet.nodes.size() << "]=";
    //            for (auto m = iSet.nodes.begin(); m != iSet.nodes.end(); ++m) {
    //                Log(LogSwitch::Szx::Model) << " " << *m;
    //            }
    //            Log(LogSwitch::Szx::Model) << endl;
    //        }
    //    }
    //}

    return true;
}

bool Solver::optimizeCoveringRelaxedBoolDecisionModel(Solution &sln, Arr<double> &coverWeights, Arr<tsm::Weight> &iSetWeights) {
    detectIndependentSetsInInvSubGraphs();

    ID nodeNum = input.graph().nodenum();

    auto &nodeColors(*sln.mutable_nodecolors());
    nodeColors.Clear();
    nodeColors.Resize(nodeNum, Problem::InvalidId);

    Arr<List<ID>> nodeInSets(nodeNum); // nodeInSets[n] is a list of independent sets containing node n.
    ID iSetId = 0;
    aux.independentSets.forEach([&](const tsm::Clique &iSet) {
        for (auto n = iSet.nodes.begin(); n != iSet.nodes.end(); ++n) {
            nodeInSets[*n].push_back(iSetId);
        }
        ++iSetId;
        return false;
    });

    MpSolver mp;

    // add decision variables.
    Arr<MpSolver::DecisionVar> isPicked(aux.independentSets.dataNum());
    for (ID s = 0; s < aux.independentSets.dataNum(); ++s) {
        isPicked[s] = mp.addVar(MpSolver::VariableType::Bool, 0, 1, 0);
    }
    Arr<MpSolver::DecisionVar> isCovered(nodeNum);
    for (auto n = 0; n < nodeNum; ++n) {
        isCovered[n] = mp.addVar(MpSolver::VariableType::Real, 0, 1, coverWeights[n]);
    }

    // add constraints.
    // node coverage.
    for (ID n = 0; n < nodeNum; ++n) {
        MpSolver::LinearExpr coveringSetNum;
        for (auto s = nodeInSets[n].begin(); s != nodeInSets[n].end(); ++s) {
            coveringSetNum += isPicked[*s];
        }
        mp.addConstraint(coveringSetNum >= isCovered[n]);
    }

    // set number.
    MpSolver::LinearExpr pickedSetNum;
    for (ID s = 0; s < aux.independentSets.dataNum(); ++s) {
        pickedSetNum += isPicked[s];
    }
    //mp.addConstraint(pickedSetNum == input.colornum());
    mp.addConstraint(pickedSetNum <= input.colornum());

    // set objective.
    mp.setOptimaOrientation(MpSolver::OptimaOrientation::Maximize);

    // solve model.
    mp.setOutput(true);
    //mp.setMaxThread(1);
    mp.setTimeLimitInSecond(1800);
    //mp.setMipFocus(MpSolver::MipFocusMode::ImproveFeasibleSolution);

    // record decision.
    if (mp.optimize()) {
        ID coloredNodeNum = 0;
        ID color = 0;
        ID s = 0;
        aux.independentSets.forEach([&](const tsm::Clique &iSet) {
            if (!mp.isTrue(isPicked[s++])) { return false; }
            for (auto n = iSet.nodes.begin(); n != iSet.nodes.end(); ++n) {
                if (nodeColors[*n] >= 0) { continue; }
                nodeColors[*n] = color;
                ++coloredNodeNum;
            }
            return (++color >= input.colornum()); // stop assigning colors if the colors are used up.
        });
        if (coloredNodeNum >= nodeNum) { return true; }
        Log(LogSwitch::Szx::Model) << "Uncovered[" << (nodeNum - coloredNodeNum) << "]=";
    }

    return false;
}

bool Solver::optimizeLocalSearch(Solution &sln) {
    ID nodeNum = input.graph().nodenum();
    ID colorNum = input.colornum();

    // cache for neighborhood moves.
    struct Move {
        ID node;
        ID color;
        //ID delta;
    };
    PriorityQueue<Move> moveQueue(2 * nodeNum, nodeNum); // moveQueue.top() is the best imrpovement neighborhood move.
    Arr2D<ID> adjColorNums(nodeNum, colorNum); // adjColorNums[n][c] is the number of adjacent nodes of node n in color c.

    // solution representation.
    struct Coloring {
        Coloring(ID nodeNumber, ID conflictNumber) : conflictNum(conflictNumber), colors(nodeNumber) {}
        Coloring(ID nodeNumber) : Coloring(nodeNumber, nodeNumber * nodeNumber) {}
        ID conflictNum;
        List<ID> colors;
    };
    Coloring optSln(nodeNum); // best solution.
    Coloring localOptSln(nodeNum); // best solution in current trajectory.
    Coloring curSln(nodeNum, 0); // current solution.
    ID &conflictNum(curSln.conflictNum);
    List<ID> &colors(curSln.colors);

    auto retreiveSln = [&](const Coloring &solution) {
        for (auto n = solution.colors.begin(); n != solution.colors.end(); ++n) { sln.add_nodecolors(*n); }
        return true;
    };

    // perturbation.
    List<ID> conflictNodes;
    conflictNodes.reserve(nodeNum);

    // reduction data.
    ID fixedNodeNum = 0;
    Arr<bool> isFixed(nodeNum, false);

    // random initialization.
    for (auto n = colors.begin(); n != colors.end(); ++n) { *n = rand.pick(colorNum); }
    for (auto n = aux.clique.nodes.begin(); n != aux.clique.nodes.end(); ++n) {
        isFixed[*n] = true;
        colors[*n] = fixedNodeNum++;
    }
    // init cache and objective.
    auto initCacheAndObj = [&]() {
        conflictNum = 0;
        moveQueue.clear();
        adjColorNums.reset(Arr2D<ID>::ResetOption::AllBits0);
        for (ID n = 0; n < nodeNum; ++n) {
            const auto &adjColorNum(adjColorNums[n]);
            if (isFixed[n]) {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    if (colors[*m] == colors[n]) { ++conflictNum; }
                }
            } else {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    ++adjColorNum[colors[*m]];
                }
                ID curConflict = adjColorNum[colors[n]];
                if (curConflict <= 0) { continue; }
                for (ID c = 0; c < colorNum; ++c) {
                    if (c == colors[n]) { continue; }
                    ID delta = adjColorNum[c] - curConflict;
                    moveQueue.push({ n, c }, delta);
                }
                conflictNum += curConflict;
            }
        }
        conflictNum /= 2;
    };

    const Iteration MaxStagnation = static_cast<Iteration>(Configuration::MaxStagnationCoefOnNodeNum * nodeNum);
    const Iteration MaxPerturbation = static_cast<Iteration>(Configuration::MaxPerterbationCoefOnNodeNum * nodeNum);
    const Iteration PerturbedNodeNum = static_cast<Iteration>(Configuration::PerturbedNodeRatio * nodeNum);
    Iteration iter = 0;
    for (Iteration perturbation = MaxPerturbation; perturbation > 0; --perturbation) {
        initCacheAndObj();
        for (Iteration stagnation = MaxStagnation; !timer.isTimeOut() && (stagnation > 0); ++iter, --stagnation) {
            // find best move.
            Move move;
            ID delta;
            ID oldColor;
            for (;;) {
                if (moveQueue.empty()) { return false; }
                delta = moveQueue.topPriority();
                moveQueue.pop(move, rand);
                oldColor = colors[move.node];
                if (move.color == oldColor) { continue; }
                ID realDelta = adjColorNums.at(move.node, move.color) - adjColorNums.at(move.node, oldColor);
                if (delta == realDelta) { break; }
            }

            // update solution.
            colors[move.node] = move.color;
            conflictNum += delta;
            Log(LogSwitch::Szx::TabuSearch) << "iter=" << iter << " opt=" << optSln.conflictNum << " cur=" << conflictNum << " node=" << move.node << " color=" << move.color << " delta=" << delta << endl;
            // update optima.
            if (conflictNum < localOptSln.conflictNum) {
                if (conflictNum <= 0) { return retreiveSln(curSln); }
                if (conflictNum < optSln.conflictNum) { optSln = curSln; }
                localOptSln = curSln;
                stagnation = MaxStagnation;
                perturbation = MaxPerturbation;
            }
            // update cache.
            for (auto n = aux.adjList[move.node].begin(); n != aux.adjList[move.node].end(); ++n) {
                if (isFixed[*n]) { continue; }
                const auto &adjColorNum(adjColorNums[*n]);
                --adjColorNum[oldColor];
                ++adjColorNum[move.color];
                ID color = colors[*n];
                ID curConflict = adjColorNum[color];
                if (curConflict <= 0) { continue; }
                if (color != oldColor) { moveQueue.push({ *n, oldColor }, adjColorNum[oldColor] - curConflict); }
                if (color != move.color) { moveQueue.push({ *n, move.color }, adjColorNum[move.color] - curConflict); }
            }
            moveQueue.push({ move.node, oldColor }, -delta);
        }

        // pertrubation.
        curSln = (rand.isPicked(1, 4) ? optSln : localOptSln); // TODO[szx][5]: parameterize the constant!
        conflictNodes.clear();
        for (ID n = 0; n < nodeNum; ++n) {
            if (adjColorNums.at(n, colors[n]) > 0) { conflictNodes.push_back(n); }
        }
        ID perturbOnConflictNodes = min(static_cast<ID>(conflictNodes.size()), static_cast<ID>(PerturbedNodeNum * Configuration::PerturbedConflictNodeRatio));
        for (ID i = 0; i < perturbOnConflictNodes; ++i) {
            ID node = conflictNodes[rand.pick(static_cast<ID>(conflictNodes.size()))];
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
        for (ID i = perturbOnConflictNodes; i < PerturbedNodeNum; ++i) {
            ID node = rand.pick(nodeNum);
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
    }

    return retreiveSln(optSln);
}

bool Solver::optimizeTabuSearchPQ(Solution &sln) {
    ID nodeNum = input.graph().nodenum();
    ID colorNum = input.colornum();
    if (iteration <= 0) { iteration = env.maxIter; }

    // cache for neighborhood moves.
    struct Move {
        ID node;
        ID color;
        //ID delta;
    };
    struct MoveEx : public Move {
        ID delta = Problem::MaxConflictNum;
        ID oldColor;
    };
    PriorityQueue<Move> moveQueue(2 * nodeNum, nodeNum); // moveQueue.top() is the best imrpovement neighborhood move.
    Arr2D<ID> adjColorNums(nodeNum, colorNum); // adjColorNums[n][c] is the number of adjacent nodes of node n in color c.

    // solution representation.
    struct Coloring {
        Coloring(ID nodeNumber, ID conflictNumber) : conflictNum(conflictNumber), colors(nodeNumber) {}
        Coloring(ID nodeNumber) : Coloring(nodeNumber, nodeNumber * nodeNumber) {}
        ID conflictNum;
        List<ID> colors;
    };
    Coloring optSln(nodeNum); // best solution.
    Coloring localOptSln(nodeNum); // best solution in current trajectory.
    Coloring curSln(nodeNum, 0); // current solution.
    ID &conflictNum(curSln.conflictNum);
    List<ID> &colors(curSln.colors);

    // tabu.
    const Iteration MaxTabuTenure = nodeNum;
    LoopQueue<List<Move>> tabuMoves(MaxTabuTenure);
    PriorityQueue<Move> tabuMoveQueue(2 * nodeNum, nodeNum);
    Arr2D<Iteration> tabuTab(nodeNum, colorNum);

    // perturbation.
    List<ID> conflictNodes;
    conflictNodes.reserve(nodeNum);

    // reduction data.
    ID fixedNodeNum = 0;
    Arr<bool> isFixed(nodeNum, false);

    // random initialization.
    for (auto n = colors.begin(); n != colors.end(); ++n) { *n = rand.pick(colorNum); }
    for (auto n = aux.clique.nodes.begin(); n != aux.clique.nodes.end(); ++n) {
        isFixed[*n] = true;
        colors[*n] = fixedNodeNum++;
    }
    // init cache and objective.
    auto initCacheAndObj = [&]() {
        conflictNum = 0;
        moveQueue.clear();
        adjColorNums.reset(Arr2D<ID>::ResetOption::AllBits0);
        tabuTab.reset(Arr2D<Iteration>::ResetOption::AllBits0);
        for (ID n = 0; n < nodeNum; ++n) {
            const auto &adjColorNum(adjColorNums[n]);
            if (isFixed[n]) {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    if (colors[*m] == colors[n]) { ++conflictNum; }
                }
            } else {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    ++adjColorNum[colors[*m]];
                }
                ID curConflict = adjColorNum[colors[n]];
                if (curConflict <= 0) { continue; }
                for (ID c = 0; c < colorNum; ++c) {
                    if (c == colors[n]) { continue; }
                    ID delta = adjColorNum[c] - curConflict;
                    moveQueue.push({ n, c }, delta);
                }
                conflictNum += curConflict;
            }
        }
        conflictNum /= 2;
    };
    
    const Iteration MaxStagnation = static_cast<Iteration>(Configuration::MaxStagnationCoefOnNodeNum * nodeNum);
    const Iteration MaxPerturbation = static_cast<Iteration>(Configuration::MaxPerterbationCoefOnNodeNum * nodeNum);
    const Iteration PerturbedNodeNum = static_cast<Iteration>(Configuration::PerturbedNodeRatio * nodeNum);
    Iteration iter = 0;
    auto isTabued = [&](const Move &move) { return iter < tabuTab.at(move.node, move.color); };
    auto pushMove = [&](const Move &move, ID objDelta) { (isTabued(move) ? tabuMoveQueue : moveQueue).push(move, objDelta); };
    Timer::TimePoint searchBegin = Timer::Clock::now();
    auto retreiveSln = [&](const Coloring &solution) {
        Math::updateMin(iteration, iter);
        double seconds = Timer::durationInSecond(searchBegin, Timer::Clock::now());
        Log(LogSwitch::Szx::TabuSearch) << "speed=" << iteration << "/" << seconds << "=" << (iteration / seconds) << endl;
        for (auto n = solution.colors.begin(); n != solution.colors.end(); ++n) { sln.add_nodecolors(*n); }
        return true;
    };
    for (Iteration perturbation = MaxPerturbation; perturbation > 0; --perturbation) {
        initCacheAndObj();
        for (Iteration stagnation = MaxStagnation; !timer.isTimeOut() && (stagnation > 0); ++iter, --stagnation) {
            // find best move.
            MoveEx move;
            MoveEx tabuMove;
            while (!moveQueue.empty()) { // best non-tabu move.
                move.delta = moveQueue.topPriority();
                moveQueue.pop(move, rand);
                move.oldColor = colors[move.node];
                if (move.color == move.oldColor) { continue; }
                if (adjColorNums.at(move.node, move.oldColor) <= 0) { continue; }
                if (isTabued(move)) { continue; }
                ID realDelta = adjColorNums.at(move.node, move.color) - adjColorNums.at(move.node, move.oldColor);
                if (move.delta == realDelta) { break; }
            }
            while (!tabuMoveQueue.empty()) { // best tabu move.
                tabuMove.delta = tabuMoveQueue.topPriority();
                if ((move.delta < Problem::MaxConflictNum)
                    && (((conflictNum + tabuMove.delta) >= localOptSln.conflictNum) || (tabuMove.delta >= move.delta))) { break; }
                tabuMoveQueue.pop(tabuMove, rand);
                tabuMove.oldColor = colors[tabuMove.node];
                if (tabuMove.color == tabuMove.oldColor) { continue; }
                if (adjColorNums.at(tabuMove.node, tabuMove.oldColor) <= 0) { continue; }
                ID realDelta = adjColorNums.at(tabuMove.node, tabuMove.color) - adjColorNums.at(tabuMove.node, tabuMove.oldColor);
                if (tabuMove.delta == realDelta) { move = tabuMove; break; }
            }

            // update solution (assuming there is always a valid move).
            colors[move.node] = move.color;
            conflictNum += move.delta;
            // update optima.
            if (conflictNum < localOptSln.conflictNum) {
                Log(LogSwitch::Szx::TabuSearch) << "iter=" << iter << " opt=" << optSln.conflictNum << " cur=" << conflictNum << endl;
                if (conflictNum <= 0) { return retreiveSln(curSln); }
                if (conflictNum < optSln.conflictNum) { optSln = curSln; }
                localOptSln = curSln;
                stagnation = MaxStagnation;
                perturbation = MaxPerturbation;
            }
            // update cache.
            for (auto n = aux.adjList[move.node].begin(); n != aux.adjList[move.node].end(); ++n) {
                if (isFixed[*n]) { continue; }
                const auto &adjColorNum(adjColorNums[*n]);
                --adjColorNum[move.oldColor];
                ++adjColorNum[move.color];
                ID color = colors[*n];
                ID curConflict = adjColorNum[color];
                if (curConflict <= 0) { continue; }
                if ((color != move.oldColor) && (color != move.color)) {
                    pushMove({ *n, move.oldColor }, adjColorNum[move.oldColor] - curConflict);
                    pushMove({ *n, move.color }, adjColorNum[move.color] - curConflict);
                } else { // the current conflict of this adjacent node is changed, so refresh the delta for moving to every other color.
                    for (ID c = 0; c < colorNum; ++c) {
                        if (c != color) { pushMove({ *n, c }, adjColorNum[c] - curConflict); }
                    }
                }
            }
            const auto &adjColorNum(adjColorNums[move.node]);
            ID curConflict = adjColorNum[move.color];
            for (ID c = 0; c < colorNum; ++c) { // the current conflict of move.node is changed, so refresh the delta for moving to every other color.
                if ((c == move.color) || (c == move.oldColor)) { continue; }
                pushMove({ move.node, c }, adjColorNum[c] - curConflict);
            }
            // update tabu.
            Move reverseMove = { move.node, move.oldColor };
            tabuMoveQueue.push(reverseMove, -move.delta);
            Iteration tabuTenure = min(conflictNum, 20) + rand.pick(2, 10); // TODO[szx][5]: parameterize the constant!
            tabuMoves.front(tabuTenure).push_back(reverseMove);
            // transfer tabu moves to non-tabu moves.
            for (auto tm = tabuMoves.front().begin(); tm != tabuMoves.front().end(); ++tm) {
                moveQueue.push(*tm, adjColorNums.at(tm->node, tm->color) - adjColorNums.at(tm->node, colors[tm->node]));
            }
            tabuMoves.front().clear();
            tabuMoves.popFront();
            tabuTab.at(move.node, move.oldColor) = iter + tabuTenure;
            //Log(LogSwitch::Szx::TabuSearch) << "iter=" << iter << " opt=" << optSln.conflictNum << " cur=" << conflictNum << " node=" << move.node << " color=" << move.color << "<=" << move.oldColor << " delta=" << move.delta << " tenure=" << tabuTenure << endl;
        }

        // pertrubation.
        curSln = (rand.isPicked(1, 4) ? optSln : localOptSln); // TODO[szx][5]: parameterize the constant!
        conflictNodes.clear();
        for (ID n = 0; n < nodeNum; ++n) {
            if (adjColorNums.at(n, colors[n]) > 0) { conflictNodes.push_back(n); }
        }
        ID perturbOnConflictNodes = min(static_cast<ID>(conflictNodes.size()), static_cast<ID>(PerturbedNodeNum * Configuration::PerturbedConflictNodeRatio));
        for (ID i = 0; i < perturbOnConflictNodes; ++i) {
            ID node = conflictNodes[rand.pick(static_cast<ID>(conflictNodes.size()))];
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
        for (ID i = perturbOnConflictNodes; i < PerturbedNodeNum; ++i) {
            ID node = rand.pick(nodeNum);
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
    }

    return retreiveSln(optSln);
}

bool Solver::optimizeTabuSearch(Solution &sln) {
    ID nodeNum = input.graph().nodenum();
    ID colorNum = input.colornum();
    if (iteration <= 0) { iteration = env.maxIter; }

    // cache for neighborhood moves.
    struct Move {
        ID node;
        ID color;
        //ID delta;
    };
    struct MoveEx : public Move {
        MoveEx() {}
        MoveEx(ID nodeId, ID colorId, ID objDelta, ID oldColorId) : Move({ nodeId, colorId }), delta(objDelta), oldColor(oldColorId) {}
        ID delta = Problem::MaxConflictNum;
        ID oldColor;
    };
    Arr2D<ID> adjColorNums(nodeNum, colorNum); // adjColorNums[n][c] is the number of adjacent nodes of node n in color c.

    // solution representation.
    struct Coloring {
        Coloring(ID nodeNumber, ID conflictNumber) : conflictNum(conflictNumber), colors(nodeNumber) {}
        Coloring(ID nodeNumber) : Coloring(nodeNumber, nodeNumber * nodeNumber) {}
        ID conflictNum;
        List<ID> colors;
    };
    Coloring optSln(nodeNum); // best solution.
    Coloring localOptSln(nodeNum); // best solution in current trajectory.
    Coloring curSln(nodeNum, 0); // current solution.
    ID &conflictNum(curSln.conflictNum);
    List<ID> &colors(curSln.colors);

    // tabu.
    const Iteration MaxTabuTenure = nodeNum;
    Arr2D<Iteration> tabuTab(nodeNum, colorNum);

    // perturbation.
    List<ID> conflictNodes;
    conflictNodes.reserve(nodeNum);

    // reduction data.
    ID fixedNodeNum = 0;
    Arr<bool> isFixed(nodeNum, false);

    // random initialization.
    for (auto n = colors.begin(); n != colors.end(); ++n) { *n = rand.pick(colorNum); }
    for (auto n = aux.clique.nodes.begin(); n != aux.clique.nodes.end(); ++n) {
        isFixed[*n] = true;
        colors[*n] = fixedNodeNum++;
    }
    // init cache and objective.
    auto initCacheAndObj = [&]() {
        conflictNum = 0;
        adjColorNums.reset(Arr2D<ID>::ResetOption::AllBits0);
        tabuTab.reset(Arr2D<Iteration>::ResetOption::AllBits0);
        for (ID n = 0; n < nodeNum; ++n) {
            const auto &adjColorNum(adjColorNums[n]);
            if (isFixed[n]) {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    if (colors[*m] == colors[n]) { ++conflictNum; }
                }
            } else {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    ++adjColorNum[colors[*m]];
                }
                conflictNum += adjColorNum[colors[n]];
            }
        }
        conflictNum /= 2;
    };

    const Iteration MaxStagnation = static_cast<Iteration>(Configuration::MaxStagnationCoefOnNodeNum * nodeNum);
    const Iteration MaxPerturbation = static_cast<Iteration>(Configuration::MaxPerterbationCoefOnNodeNum * nodeNum);
    const Iteration PerturbedNodeNum = static_cast<Iteration>(Configuration::PerturbedNodeRatio * nodeNum);
    Iteration iter = 0;
    Timer::TimePoint searchBegin = Timer::Clock::now();
    auto retreiveSln = [&](const Coloring &solution) {
        Math::updateMin(iteration, iter);
        double seconds = Timer::durationInSecond(searchBegin, Timer::Clock::now());
        Log(LogSwitch::Szx::TabuSearch) << "speed=" << iteration << "/" << seconds << "=" << (iteration / seconds) << endl;
        for (auto n = solution.colors.begin(); n != solution.colors.end(); ++n) { sln.add_nodecolors(*n); }
        return true;
    };
    for (Iteration perturbation = MaxPerturbation; perturbation > 0; --perturbation) {
        initCacheAndObj();
        for (Iteration stagnation = MaxStagnation; !timer.isTimeOut() && (stagnation > 0); ++iter, --stagnation) {
            // find best move.
            MoveEx move;
            MoveEx tabuMove;
            Sampling1 picker(rand, Sampling1::StartCount::WithPresetElement);
            for (ID n = 0; n < nodeNum; ++n) {
                if (isFixed[n]) { continue; }
                const auto &adjColorNum(adjColorNums[n]);
                if (adjColorNum[colors[n]] <= 0) { continue; }
                for (ID c = 0; c < colorNum; ++c) {
                    if (c == colors[n]) { continue; }
                    ID delta = adjColorNums.at(n, c) - adjColorNums.at(n, colors[n]);
                    if (iter < tabuTab.at(n, c)) { // tabu move.
                        if (picker.isMinimal(tabuMove.delta, delta)) {
                            tabuMove = MoveEx(n, c, delta, colors[n]);
                        }
                    } else { // non-tabu move.
                        if (picker.isMinimal(move.delta, delta)) {
                            move = MoveEx(n, c, delta, colors[n]);
                        }
                    }
                }
            }
            if (((conflictNum + tabuMove.delta) < localOptSln.conflictNum) && (tabuMove.delta < move.delta)) { move = tabuMove; }
            if (move.delta >= Problem::MaxConflictNum) { move = tabuMove; }

            // update solution (assuming there is always a valid move).
            colors[move.node] = move.color;
            conflictNum += move.delta;

            // update optima.
            if (conflictNum < localOptSln.conflictNum) {
                Log(LogSwitch::Szx::TabuSearch) << "iter=" << iter << " opt=" << optSln.conflictNum << " cur=" << conflictNum << endl;
                if (conflictNum <= 0) { return retreiveSln(curSln); }
                if (conflictNum < optSln.conflictNum) { optSln = curSln; }
                localOptSln = curSln;
                stagnation = MaxStagnation;
                perturbation = MaxPerturbation;
            }
            // update cache.
            for (auto n = aux.adjList[move.node].begin(); n != aux.adjList[move.node].end(); ++n) {
                if (isFixed[*n]) { continue; }
                const auto &adjColorNum(adjColorNums[*n]);
                --adjColorNum[move.oldColor];
                ++adjColorNum[move.color];
            }
            // update tabu.
            Iteration tabuTenure = min(conflictNum, 20) + rand.pick(2, 10); // TODO[szx][5]: parameterize the constant!
            tabuTab.at(move.node, move.oldColor) = iter + tabuTenure;
            //Log(LogSwitch::Szx::TabuSearch) << "iter=" << iter << " opt=" << optSln.conflictNum << " cur=" << conflictNum << " node=" << move.node << " color=" << move.color << "<=" << move.oldColor << " delta=" << move.delta << " tenure=" << tabuTenure << endl;
        }

        // pertrubation.
        curSln = (rand.isPicked(1, 4) ? optSln : localOptSln); // TODO[szx][5]: parameterize the constant!
        conflictNodes.clear();
        for (ID n = 0; n < nodeNum; ++n) {
            if (adjColorNums.at(n, colors[n]) > 0) { conflictNodes.push_back(n); }
        }
        ID perturbOnConflictNodes = min(static_cast<ID>(conflictNodes.size()), static_cast<ID>(PerturbedNodeNum * Configuration::PerturbedConflictNodeRatio));
        for (ID i = 0; i < perturbOnConflictNodes; ++i) {
            ID node = conflictNodes[rand.pick(static_cast<ID>(conflictNodes.size()))];
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
        for (ID i = perturbOnConflictNodes; i < PerturbedNodeNum; ++i) {
            ID node = rand.pick(nodeNum);
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
    }

    return retreiveSln(optSln);
}

bool Solver::optimizeAuctionSearch(Solution &sln) {
    ID nodeNum = input.graph().nodenum();
    ID colorNum = input.colornum();

    // cache for neighborhood moves.
    struct Move {
        ID node;
        ID color;
        //ID delta;
    };
    // OPTIMIZE[szx][0]: use L2BucketQueue?
    PriorityQueue<Move> moveQueue(2 * nodeNum * Configuration::ConflictWeightBase, nodeNum * Configuration::ConflictWeightBase); // moveQueue.top() is the best imrpovement neighborhood move.
    Arr2D<ID> adjColorNums(nodeNum, colorNum); // adjColorNums[n][c] is the number of adjacent nodes of node n in color c.

    // solution representation.
    struct Coloring {
        Coloring(ID nodeNumber, ID conflictNumber) : conflictNum(conflictNumber),
            weightedConflict(conflictNumber * Configuration::ConflictWeightBase), colors(nodeNumber) {
        }
        Coloring(ID nodeNumber) : Coloring(nodeNumber, nodeNumber * nodeNumber) {}
        ID conflictNum;
        ID weightedConflict;
        List<ID> colors;
    };
    Coloring optSln(nodeNum); // best solution.
    Coloring localOptSln(nodeNum); // best solution in current trajectory.
    Coloring curSln(nodeNum, 0); // current solution.
    ID &conflictNum(curSln.conflictNum);
    ID &weightedConflict(curSln.weightedConflict);
    List<ID> &colors(curSln.colors);

    auto retreiveSln = [&](const Coloring &solution) {
        for (auto n = solution.colors.begin(); n != solution.colors.end(); ++n) { sln.add_nodecolors(*n); }
        return true;
    };

    // weight.
    Arr2D<ID> weights(nodeNum, colorNum, Configuration::ConflictWeightBase);
    auto calcWeight = [&](ID node, ID color) {
        return weights.at(node, color) * adjColorNums.at(node, color);
    };

    // perturbation.
    List<ID> conflictNodes;
    conflictNodes.reserve(nodeNum);

    // reduction data.
    ID fixedNodeNum = 0;
    Arr<bool> isFixed(nodeNum, false);

    // random initialization.
    for (auto n = colors.begin(); n != colors.end(); ++n) { *n = rand.pick(colorNum); }
    for (auto n = aux.clique.nodes.begin(); n != aux.clique.nodes.end(); ++n) {
        isFixed[*n] = true;
        colors[*n] = fixedNodeNum++;
    }
    // init cache and objective.
    auto initCacheAndObj = [&]() {
        conflictNum = 0;
        moveQueue.clear();
        adjColorNums.reset(Arr2D<ID>::ResetOption::AllBits0);
        for (ID n = 0; n < nodeNum; ++n) {
            const auto &adjColorNum(adjColorNums[n]);
            if (isFixed[n]) {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    if (colors[*m] == colors[n]) { ++conflictNum; }
                }
            } else {
                for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                    ++adjColorNum[colors[*m]];
                }
                ID curConflict = adjColorNum[colors[n]]; // calcWeight(n, colors[n]);
                if (curConflict <= 0) { continue; }
                for (ID c = 0; c < colorNum; ++c) {
                    if (c == colors[n]) { continue; }
                    ID delta = Configuration::ConflictWeightBase * (adjColorNum[c] - curConflict); // calcWeight(n, c);
                    moveQueue.push({ n, c }, delta);
                }
                conflictNum += curConflict;
            }
        }
        conflictNum /= 2;
        weightedConflict = conflictNum * Configuration::ConflictWeightBase;
    };

    const Iteration MaxStagnation = static_cast<Iteration>(Configuration::MaxStagnationCoefOnNodeNum * nodeNum);
    const Iteration MaxPerturbation = static_cast<Iteration>(Configuration::MaxPerterbationCoefOnNodeNum * nodeNum);
    const Iteration PerturbedNodeNum = static_cast<Iteration>(Configuration::PerturbedNodeRatio * nodeNum);
    Iteration iter = 0;
    for (Iteration perturbation = MaxPerturbation; perturbation > 0; --perturbation) {
        initCacheAndObj();
        for (Iteration stagnation = MaxStagnation; !timer.isTimeOut() && (stagnation > 0); ++iter, --stagnation) {
            // find best move.
            Move move;
            ID delta;
            ID oldColor;
            for (;;) {
                if (moveQueue.empty()) { return false; }
                delta = moveQueue.topPriority();
                moveQueue.pop(move, rand);
                oldColor = colors[move.node];
                if (move.color == oldColor) { continue; }
                ID realDelta = calcWeight(move.node, move.color) - calcWeight(move.node, oldColor);
                if (delta == realDelta) { break; }
            }

            // update solution.
            colors[move.node] = move.color;
            weightedConflict += delta;
            (conflictNum += adjColorNums.at(move.node, move.color)) -= adjColorNums.at(move.node, oldColor);
            Log(LogSwitch::Szx::TabuSearch) << "iter=" << iter << "opt=" << optSln.conflictNum << " cur=" << conflictNum << " w=" << weightedConflict << " node=" << move.node << " color=" << move.color << " delta=" << delta << endl;
            // update optima.
            if (conflictNum < localOptSln.conflictNum) {
                if (conflictNum <= 0) { return retreiveSln(curSln); }
                if (conflictNum < optSln.conflictNum) { optSln = curSln; }
                localOptSln = curSln;
                stagnation = MaxStagnation;
                perturbation = MaxPerturbation;
            }
            // update weights when trapped in local optima.
            if (delta >= 0) {
                ++weights.at(move.node, oldColor);
                //weightedConflict += adjColorNums.at(move.node, oldColor);
            }
            // update cache.
            for (auto n = aux.adjList[move.node].begin(); n != aux.adjList[move.node].end(); ++n) {
                if (isFixed[*n]) { continue; }
                --adjColorNums.at(*n, oldColor);
                ++adjColorNums.at(*n, move.color);
                ID color = colors[*n];
                ID curWeightedConflict = calcWeight(*n, color);
                if (curWeightedConflict <= 0) { continue; }
                if (color != oldColor) { moveQueue.push({ *n, oldColor }, calcWeight(*n, oldColor) - curWeightedConflict); }
                if (color != move.color) { moveQueue.push({ *n, move.color }, calcWeight(*n, move.color) - curWeightedConflict); }
            }
            moveQueue.push({ move.node, oldColor }, -delta);
        }

        // pertrubation.
        curSln = (rand.isPicked(1, 4) ? optSln : localOptSln); // TODO[szx][5]: parameterize the constant!
        conflictNodes.clear();
        for (ID n = 0; n < nodeNum; ++n) {
            if (adjColorNums.at(n, colors[n]) > 0) { conflictNodes.push_back(n); }
        }
        ID perturbOnConflictNodes = min(static_cast<ID>(conflictNodes.size()), static_cast<ID>(PerturbedNodeNum * Configuration::PerturbedConflictNodeRatio));
        for (ID i = 0; i < perturbOnConflictNodes; ++i) {
            ID node = conflictNodes[rand.pick(static_cast<ID>(conflictNodes.size()))];
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
        for (ID i = perturbOnConflictNodes; i < PerturbedNodeNum; ++i) {
            ID node = rand.pick(nodeNum);
            if (!isFixed[node]) { colors[node] = rand.pick(colorNum); }
        }
    }

    return retreiveSln(optSln);
}

bool Solver::optimizeCliqueReduction(Solution &sln) {
    auto &nodeColors(*sln.mutable_nodecolors());
    nodeColors.Resize(input.graph().nodenum(), Problem::InvalidId);

    // node mapping.
    ID colorNum = input.colornum();
    ID nodeNum = input.graph().nodenum();
    ID virtualNodeNum = 0;
    Arr2D<ID> indices(nodeNum, colorNum); // virtual node id if it exists.
    indices.reset(Arr2D<ID>::ResetOption::AllBits1);
    for (ID n = 0; n < nodeNum; ++n) {
        if (aux.fixedColors[n] > Problem::InvalidId) { continue; } // eliminate nodes with fixed colors.
        List<bool> skipColors(colorNum, false); // eliminate colors which is used by adjacent nodes with fixed colors.
        for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
            if (aux.fixedColors[*m] > Problem::InvalidId) { skipColors[aux.fixedColors[*m]] = true; }
        }
        for (ID c = 0; c < colorNum; ++c) {
            if (!skipColors[c]) { indices.at(n, c) = virtualNodeNum++; }
        }
    }
    Log(LogSwitch::CliqueReduction) << "RemainingNodeNum=" << virtualNodeNum << "/" << nodeNum * colorNum << endl;

    // generate virtual graph.
    tsm::AdjMat adjMat(virtualNodeNum, virtualNodeNum);
    adjMat.reset(tsm::AdjMat::ResetOption::AllBits0);

    ID maxEdgeNum = virtualNodeNum * (virtualNodeNum - 1) / 2;
    ID edgeNum = maxEdgeNum;
    for (ID n = 0; n < nodeNum; ++n) {
        if (aux.fixedColors[n] > Problem::InvalidId) { continue; }
        for (ID c = 0; c < colorNum; ++c) {
            if (indices.at(n, c) <= Problem::InvalidId) { continue; }
            // single color.
            for (ID k = 0; k < c; ++k) {
                if (indices.at(n, k) <= Problem::InvalidId) { continue; }
                adjMat.at(indices.at(n, c), indices.at(n, k)) = true;
                adjMat.at(indices.at(n, k), indices.at(n, c)) = true;
                --edgeNum;
            }
            // conflict avoidance.
            for (auto m = aux.adjList[n].begin(); m != aux.adjList[n].end(); ++m) {
                if (aux.fixedColors[*m] > Problem::InvalidId) { continue; }
                if (indices.at(*m, c) <= Problem::InvalidId) { continue; }
                adjMat.at(indices.at(n, c), indices.at(*m, c)) = true;
                adjMat.at(indices.at(*m, c), indices.at(n, c)) = true;
                --edgeNum;
            }
        }
    }
    Log(LogSwitch::CliqueReduction) << "VirtualEdgeNum=" << edgeNum << "/" << maxEdgeNum << endl;

    // solve.
    tsm::Clique iSet;
    tsm::solveWeightedIndependentSet(iSet, adjMat, tsm::Forever, nodeNum - static_cast<ID>(aux.clique.nodes.size()));

    // retrive solution.
    sort(iSet.nodes.begin(), iSet.nodes.end());
    auto vn = iSet.nodes.begin();
    for (ID n = 0; n < nodeNum; ++n) {
        if (aux.fixedColors[n] > Problem::InvalidId) {
            nodeColors[n] = aux.fixedColors[n];
        } else {
            for (ID c = 0; c < colorNum; ++c) {
                if (*vn == indices.at(n, c)) { nodeColors[n] = c; break; }
            }
            ++vn;
        }
    }

    return false;
}

void Solver::detectClique() {
    if (!aux.fixedColors.empty()) { return; } // already initialized.

    aux.fixedColors.resize(input.graph().nodenum(), Problem::InvalidId);

    // OPTIMIZE[szx][5]: the tsm seems to find the max clique with max degree automatically?
    Arr<tsm::Weight> weights(input.graph().nodenum()); // a node with greater degree is more important.
    for (ID n = 0; n < input.graph().nodenum(); ++n) {
        weights[n] = 1; // Configuration::TsmPrecision + static_cast<tsm::Weight>(aux.adjList[n].size());
    }

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
    }, weights, cfg.msCliqueDetectionTimeout);

    Log(LogSwitch::Szx::Preprocess) << "Clique[" << aux.clique.weight << "]=";
    ID c = 0;
    for (auto n = aux.clique.nodes.begin(); n != aux.clique.nodes.end(); ++n) {
        aux.fixedColors[*n] = c++;
        Log(LogSwitch::Szx::Preprocess) << " " << *n;
    }
    Log(LogSwitch::Szx::Preprocess) << endl;
}

void Solver::detectIndependentSet() {
    if (aux.independentSets.dataNum() > 0) { return; } // already initialized.

    ID nodeNum = input.graph().nodenum();

    aux.independentSets = CombinationMap<tsm::Clique>(nodeNum);

    tsm::Clique iSet;
    tsm::solveWeightedIndependentSet(iSet, aux.adjMat, cfg.msCliqueDetectionTimeout);
    sort(iSet.nodes.begin(), iSet.nodes.end());
    aux.independentSets.set(iSet.nodes, iSet);
    Log(LogSwitch::Szx::Preprocess) << "IndependentSet[" << iSet.weight << "]=";
    for (auto n = iSet.nodes.begin(); n != iSet.nodes.end(); ++n) {
        Log(LogSwitch::Szx::Preprocess) << " " << *n;
    }
    Log(LogSwitch::Szx::Preprocess) << endl;
}

void Solver::detectIndependentSetsInInvSubGraphs() {
    if (!aux.subInvGraphs.empty()) { return; } // already initialized.
    detectIndependentSet();

    ID nodeNum = input.graph().nodenum();

    tsm::Clique iSet;
    aux.subInvGraphs.resize(nodeNum);
    for (ID n = 0; n < nodeNum; ++n) {
        Subgraph &subInvGraph(aux.subInvGraphs[n]);
        for (ID m = 0; m < nodeNum; ++m) {
            if (!aux.adjMat.at(n, m)) { subInvGraph.idMap.push_back(m); } // including node n itself.
        }

        ID subInvGraphNodeNum = static_cast<ID>(subInvGraph.idMap.size());
        subInvGraph.adjMat.init(subInvGraphNodeNum, subInvGraphNodeNum); // no default value since every item will be set.
        for (ID i = 0; i < subInvGraphNodeNum; ++i) {
            for (ID j = 0; j < subInvGraphNodeNum; ++j) {
                //if (i == j) { continue; } // the tsm will ignore edge to itself.
                subInvGraph.adjMat.at(i, j) = !aux.adjMat.at(subInvGraph.idMap[i], subInvGraph.idMap[j]);
            }
        }
        // OPTIMIZE[szx][5]: record the combination of removed nodes to avoid computing clique on the same sub-graph.
        tsm::solveWeightedMaxClique(iSet, subInvGraph.adjMat, cfg.msCliqueDetectionTimeout);
        subInvGraph.mapBack(iSet.nodes);
        sort(iSet.nodes.begin(), iSet.nodes.end());
        if (aux.independentSets.set(iSet.nodes, iSet)) {
            Log(LogSwitch::Szx::Preprocess) << n << ".IndependentSet[" << iSet.weight << "]=";
            for (auto m = iSet.nodes.begin(); m != iSet.nodes.end(); ++m) {
                Log(LogSwitch::Szx::Preprocess) << " " << *m;
            }
            Log(LogSwitch::Szx::Preprocess) << endl;
        }
    }
}

void Solver::detectIndependentSetsInSubGraphs() {
    if (!aux.subGraphs.empty()) { return; } // already initialized.

    ID nodeNum = input.graph().nodenum();

    tsm::Clique iSet;
    aux.subGraphs.resize(nodeNum);
    aux.independentSetInSubGraphs.resize(nodeNum);
    for (ID n = 0; n < nodeNum; ++n) {
        Subgraph &subGraph(aux.subGraphs[n]);
        for (ID m = 0; m < nodeNum; ++m) {
            if (aux.adjMat.at(n, m)) { subGraph.idMap.push_back(m); } // not including node n itself.
        }

        ID subGraphNodeNum = static_cast<ID>(subGraph.idMap.size());
        subGraph.adjMat.init(subGraphNodeNum, subGraphNodeNum); // no default value since every item will be set.
        for (ID i = 0; i < subGraphNodeNum; ++i) {
            for (ID j = 0; j < subGraphNodeNum; ++j) {
                //if (i == j) { continue; } // the tsm will ignore edge to itself.
                subGraph.adjMat.at(i, j) = aux.adjMat.at(subGraph.idMap[i], subGraph.idMap[j]);
            }
        }
        tsm::solveWeightedIndependentSet(iSet, subGraph.adjMat, cfg.msCliqueDetectionTimeout);
        subGraph.mapBack(iSet.nodes);
        sort(iSet.nodes.begin(), iSet.nodes.end());
        aux.independentSetInSubGraphs[n] = iSet;

        Log(LogSwitch::Szx::Preprocess) << n << ".IndependentSet[" << iSet.weight << "]=";
        for (auto m = iSet.nodes.begin(); m != iSet.nodes.end(); ++m) {
            Log(LogSwitch::Szx::Preprocess) << " " << *m;
        }
        Log(LogSwitch::Szx::Preprocess) << endl;
    }
}
#pragma endregion Solver

}
