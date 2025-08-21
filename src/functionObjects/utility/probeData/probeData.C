/*---------------------------------------------------------------------------*\

                        Open Source CFD Group

\*---------------------------------------------------------------------------*/

#include "signalProcess.H"
#include "Field.H"
#include "SubField.H"
#include <fstream>
#include <math.h>
#include <limits>

namespace Foam {

Field<scalarField> readScalarProbe(fileName file, scalar begin_time, scalar& delta_t, label& count) {
    scalarField time;
    Field<scalarField> result;

    std::ifstream ifs(file.c_str());
    if (!ifs.is_open()) {
        FatalError << "Unable to open file" << nl;
        delta_t = 0;
        count = 0;
        return  {};
    } else Info << file << " opened" << nl;

    std::string line;
    label n_probes = 0;
    label line_index = 0;
    label passed_data = 0;
    bool probe_data_read = false;
    while (std::getline(ifs, line)) {
        if (line[0] == '#' && !probe_data_read) {
            n_probes++;
            continue;
        } else if (line[0] == '#' && probe_data_read) {
            Info << "Unable to parse probe data";
            return { };
        } else if (!probe_data_read) {
            n_probes -= 2;
            probe_data_read = true;
            result.resize(n_probes);
            
            auto state_backup = ifs.rdstate();
            auto pos_backup = ifs.tellg();
            label  n_lines = std::count(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>(), '\n') + 1;
            ifs.clear();
            ifs.seekg(pos_backup);
            ifs.setstate(state_backup);
            
            for (auto& r : result) r.resize(n_lines);
            time.resize(n_lines);
        }
        
        std::stringstream liness(line);

        double data;
        liness >> data;
        if (data >= begin_time) {
            time[line_index] = data;
            for (label i = 0; i < n_probes; ++i) {
                liness >> data;
                result[i][line_index] = data;
            }
            line_index++;
        } else passed_data++;
    }

    for (auto& r : result) r.resize(r.size() - passed_data);
    time.resize(time.size() - passed_data);

    if (time.size() > 1) {
        delta_t = time[1] - time[0];
        bool changing = false;
        double min_delta_t = delta_t;
        double max_delta_t = delta_t;
        double t = 0;
        double pt = 0;

        for (label i = 2; i < time.size(); ++i) {
            double dt = time[i] - time[i - 1];
            if (std::abs(delta_t - dt) > std::numeric_limits<double>::epsilon()) {
                changing = true;
                min_delta_t = std::min(min_delta_t, dt);
                max_delta_t = std::max(max_delta_t, dt);
                t = time[i];
                pt = time[i - 1];
            }
        }

        if (changing) {
            Info << "Changing timestep value found at time " << t << ", previous time is " << pt << ", ref. Delta T is " << delta_t << nl
                << "Min. Delta T = " << min_delta_t << nl
                << "Max. Delta T = " << max_delta_t << nl;
        } else {
            Info << "Data time range is [" << time[0] << ", " << time[time.size() - 1] << "]. Timesteping is " << delta_t << " with " << time.size() << " datas" << nl;
        }

        count = time.size();

        return result;
    }
    
    delta_t = 0;
    count = 0;
    return  {};
}

Field<vectorField> readVectorProbe(fileName file, scalar begin_time, scalar& delta_t, label& count) {
    scalarField time;
    Field<vectorField> result;

    std::ifstream ifs(file.c_str());
    if (!ifs.is_open()) {
        FatalError << "Unable to open file";
        delta_t = 0;
        count = 0;
        return  {};
    }

    std::string line;
    label n_probes = 0;
    label line_index = 0;
    label passed_data = 0;
    bool probe_data_read = false;
    while (std::getline(ifs, line)) {
        if (line[0] == '#' && !probe_data_read) {
            n_probes++;
            continue;
        } else if (line[0] == '#' && probe_data_read) {
            Info << "Unable to parse probe data";
            return {  };
        } else if (!probe_data_read) {
            n_probes -= 2;
            probe_data_read = true;
            result.resize(n_probes);

            auto state_backup = ifs.rdstate();
            auto pos_backup = ifs.tellg();
            label  n_lines = std::count(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>(), '\n') + 1;
            ifs.clear();
            ifs.seekg(pos_backup);
            ifs.setstate(state_backup);

            for (auto& r : result) r.resize(n_lines);
            time.resize(n_lines);
        }

        std::stringstream liness(line);

        double data;
        double u, v, w;
        liness >> data;
        if (data >= begin_time) {
            time[line_index] = data;
            for (label i = 0; i < n_probes; ++i) {
                liness.ignore(256, '(');
                std::string vector_str;
                std::getline(liness, vector_str, ')');
                std::stringstream vectorss(vector_str);
                vectorss >> u;
                vectorss >> v;
                vectorss >> w;
                result[i][line_index] = { u, v, w };
            }
            line_index++;
        } else passed_data++;
    }

    for (auto& r : result) r.resize(r.size() - passed_data);
    time.resize(time.size() - passed_data);

    if (time.size() > 1) {
        delta_t = time[1] - time[0];
        bool changing = false;
        double min_delta_t = delta_t;
        double max_delta_t = delta_t;
        double t = 0;
        double pt = 0;

        for (label i = 2; i < time.size(); ++i) {
            double dt = time[i] - time[i - 1];
            if (std::abs(delta_t - dt) > std::numeric_limits<double>::epsilon()) {
                changing = true;
                min_delta_t = std::min(min_delta_t, dt);
                max_delta_t = std::max(max_delta_t, dt);
                t = time[i];
                pt = time[i - 1];
            }
        }

        if (changing) {
            Info << "Changing timestep value found at time " << t << ", previous time is " << pt << ", ref. Delta T is " << delta_t << nl
                << "Min. Delta T = " << min_delta_t << nl
                << "Max. Delta T = " << max_delta_t << nl;
        } else {
            Info << "Data time range is [" << time[0] << ", " << time[time.size() - 1] << "]. Timesteping is " << delta_t << " with " << time.size() << " datas" << nl;
        }

        count = time.size();

        return result;
    }

    delta_t = 0;
    count = 0;
    return  {};
}

}


// ************************************************************************* //
