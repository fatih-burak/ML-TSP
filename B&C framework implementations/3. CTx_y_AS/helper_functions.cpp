#include "helper_functions.h"
constexpr double PI = 3.14159265358979323846;
constexpr double RRR = 6378.388;


double getCPUTime() {
	return (double)clock() / CLOCKS_PER_SEC;
}

void trim(std::string& s) { // to trim the strings white spaces from the left and right
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));               // ltrim
    s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);               // rtrim
}

vector<int> setDifference(const vector<int>& vertices, const set<int>& cycle) {
    /// <summary>
    /// Returns the set difference of vertices \ cycle
    /// </summary>
    vector<int> result;

    // Sort vertices (set is already sorted)
    vector<int> sorted_vertices = vertices;
    sort(sorted_vertices.begin(), sorted_vertices.end());

    // Compute set difference
    set_difference(sorted_vertices.begin(), sorted_vertices.end(),
        cycle.begin(), cycle.end(),
        back_inserter(result));

    return result;
}


double degToRad(double x) {
    int deg = (int)x;
    double min = x - deg;
    double decimal = deg + (5.0 * min) / 3.0;
    return PI * decimal / 180.0;
}

int geoDistance(double lat1, double lon1, double lat2, double lon2) {
    // Computes the great-circle (spherical) distance between two coordinates in latitude/longitude (used in GEO instances).
    double phi1 = degToRad(lat1);
    double phi2 = degToRad(lat2);
    double lambda1 = degToRad(lon1);
    double lambda2 = degToRad(lon2);

    double q1 = cos(lambda1 - lambda2);
    double q2 = cos(phi1 - phi2);
    double q3 = cos(phi1 + phi2);

    double dist = RRR * acos(0.5 * ((1 + q1) * q2 - (1 - q1) * q3)) + 1.0;
    return (int)(dist);
}

int euclideanDistance(double x1, double y1, double x2, double y2) {
    // Computes standard Euclidean distance in 2D space, rounded to the nearest integer (used in EUC_2D).
    return (int)(sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) + 0.5); // round
}

int ceilDistance(double x1, double y1, double x2, double y2) {
    // Computes 2D Euclidean distance but rounded up to the nearest integer (used in CEIL_2D).
    return (int)(ceil(sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))));
}

int attDistance(double x1, double y1, double x2, double y2) {
    // Computes pseudo-Euclidean distance used in ATT instances (approximates real-world travel costs better).
    double dx = x1 - x2;
    double dy = y1 - y2;
    double rij = sqrt((dx * dx + dy * dy) / 10.0);
    int tij = (int)(rij + 0.5);
    if (tij < rij) return tij + 1;
    else return tij;
}

int computeDistance(const string& edgeWeightType,
    double x1, double y1, double x2, double y2) {
    if (edgeWeightType == "GEO") {
        return geoDistance(x1, y1, x2, y2);
    }
    else if (edgeWeightType == "EUC_2D") {
        return euclideanDistance(x1, y1, x2, y2);
    }
    else if (edgeWeightType == "CEIL_2D") {
        return ceilDistance(x1, y1, x2, y2);
    }
    else if (edgeWeightType == "ATT") {
        return attDistance(x1, y1, x2, y2);
    }
    else {
        throw runtime_error("Unsupported EDGE_WEIGHT_TYPE: " + edgeWeightType);
    }
}


Instance readInstance(const string& filePath) {
	// define variables
	Instance instance;

    ifstream file(filePath);
    if (!file.is_open()) throw runtime_error("Unable to open file " + filePath);

    string line, section, edgeWeightType, edgeWeightFormat;
    int dimension = -1;

    // Read header
    while (getline(file, line)) {
        //cout << line << endl;
        istringstream iss(line);
        string key;
        string value;
        if (line.find(':') != string::npos) {
            getline(iss, key, ':');
            getline(iss, value);
            trim(value);
        }
        else
            getline(iss, key);

        trim(key);

        
        if (line.find("NAME") != string::npos) {
            instance.name = value;
            cout << "Intance name = " << instance.name << endl;
            continue;
        }
        else if (line.find("NODE_COORD_SECTION") != string::npos) {
            section = "NODE_COORD_SECTION";
            break;
        }
        else if (line.find("EDGE_WEIGHT_SECTION") != string::npos) {
            section = "EDGE_WEIGHT_SECTION";
            break;
        }

        if (key.find("DIMENSION") != string::npos) {
            dimension = stoi(value);
            //cout << dimension << endl;
        }
        else if (key.find("EDGE_WEIGHT_TYPE") != string::npos) {
            edgeWeightType = value;
            transform(edgeWeightType.begin(), edgeWeightType.end(), edgeWeightType.begin(), ::toupper);
            cout << edgeWeightType << endl;
        }
        else if (key.find("EDGE_WEIGHT_FORMAT") != string::npos) {
            edgeWeightFormat = value;
            //cout << edgeWeightFormat << endl;
        }

    }

    if (dimension <= 0)
        throw runtime_error("Invalid or missing DIMENSION");

    instance.nbCity = dimension;
    instance.V.resize(dimension);
    for (int i = 0; i < dimension; ++i)
        instance.V[i] = i;

    instance.dist.resize(dimension, vector<double>(dimension, 0));
    //instance.E.resize(dimension);

    if (section == "NODE_COORD_SECTION") {
        instance.xCoord.resize(dimension);
        instance.yCoord.resize(dimension);

        for (int i = 0; i < dimension; ++i) {
            int index;
            double x, y;
            file >> index >> x >> y;
            instance.xCoord[index - 1] = x;
            instance.yCoord[index - 1] = y;
        }

        // Compute distances based on EDGE_WEIGHT_TYPE
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (i == j) continue;
                instance.dist[i][j] = computeDistance(
                    edgeWeightType,
                    instance.xCoord[i], instance.yCoord[i],
                    instance.xCoord[j], instance.yCoord[j]
                );
            }
        }


    }
    else if (section == "EDGE_WEIGHT_SECTION") {

        if (edgeWeightFormat == "UPPER_ROW") {
            // Fill upper triangle
            int count = 0;
            vector<double> upper;
            double w;
            while (file >> w)
                upper.push_back(w);

            int idx = 0;
            for (int i = 0; i < dimension; ++i) {
                for (int j = i + 1; j < dimension; ++j) {
                    instance.dist[i][j] = instance.dist[j][i] = upper[idx++];
                }
            }
        }
        else if (edgeWeightFormat == "FULL_MATRIX") {
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j < dimension; ++j) {
                    double w;
                    if (!(file >> w)) {
                        throw runtime_error("Invalid FULL_MATRIX data.");
                    }
                    instance.dist[i][j] = w;
                    //if (i != j) {
                        //instance.E[i].push_back(j);
                    //}
                }
            }
        }
        else if (edgeWeightFormat == "LOWER_DIAG_ROW") {
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j <= i; ++j) {
                    double w;
                    if (!(file >> w)) {
                        throw runtime_error("Invalid LOWER_DIAG_ROW data.");
                    }
                    instance.dist[i][j] = instance.dist[j][i] = w;
                }
            }

        }
        else if (edgeWeightFormat == "UPPER_DIAG_ROW") {
            for (int i = 0; i < dimension; ++i) {
                for (int j = i; j < dimension; ++j) {
                    double w;
                    if (!(file >> w)) {
                        throw runtime_error("Invalid UPPER_DIAG_ROW data.");
                    }
                    instance.dist[i][j] = instance.dist[j][i] = w;
                }
            }

        }
        else
            throw runtime_error("Unsupported EDGE_WEIGHT_FORMAT!");

    }

    instance.startTimeSaved = getCPUTime();

    return instance;
}

void printInfo(const string& pathAndFileout, const Instance& inst, const Solution& sol) {
    string nameFile = pathAndFileout;
    std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
    cout << sol.LB << endl;
    file << fixed << setprecision(2) << inst.name << "\t" << sol.opt << "\t" << sol.timeT << "\t" << sol.iter << "\t" << sol.cuts.size() << "\t" << sol.eliminatedSECs << "\t" << sol.optimalObjVal << "\t" << sol.LB << endl;
    file.close();
}
