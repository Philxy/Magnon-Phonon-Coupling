#include "../include/couplingParameters.h"

std::vector<CouplingParameter> readCouplingParameters(const std::string &filename, Axis displacementDirection)
{
    std::vector<CouplingParameter> parameters;
    std::ifstream file(filename);
    std::string line;

    std::cout.precision(20); // does this line only determine the precision of printing or the precision of the string to double conversion?

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return parameters; // Return empty vector if file can't be opened
    }

    // Skip the first line
    std::getline(file, line);

    while (getline(file, line))
    {

        if (line.empty())
        {
            continue;
        }

        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        CouplingParameter param;
        std::istringstream iss(line);
        std::vector<double> values;
        std::string value;

        // retrieve all the values from the line and push them to 'values'
        while (std::getline(iss, value, ','))
        {
            if (!value.empty())
            { // Check if the string is not empty
                try
                {
                    values.push_back(std::stod(value)); // Convert to double and store
                }
                catch (const std::invalid_argument &e)
                {
                    std::cerr << "Invalid argument for std::stod: " << value << std::endl;
                }
                catch (const std::out_of_range &e)
                {
                    std::cerr << "Value out of range for std::stod: " << value << std::endl;
                    continue; // Handle the error
                }
            }
        }

        param.displacementDirection = displacementDirection;

        param.x_ij = values[0];
        param.y_ij = values[1];
        param.z_ij = values[2];

        param.x_jk = values[3];
        param.y_jk = values[4];
        param.z_jk = values[5];

        param.x_ki = values[6];
        param.y_ki = values[7];
        param.z_ki = values[8];

        param.J_xx = values[9];
        param.J_xy = values[10];
        param.J_xz = values[11];
        param.J_yx = values[12];
        param.J_yy = values[13];
        param.J_yz = values[14];
        param.J_zx = values[15];
        param.J_zy = values[16];
        param.J_zz = values[17];

        param.J[0][0] = values[9];
        param.J[0][1] = values[10];
        param.J[0][2] = values[11];
        param.J[1][0] = values[12];
        param.J[1][1] = values[13];
        param.J[1][2] = values[14];
        param.J[2][0] = values[15];
        param.J[2][1] = values[16];
        param.J[2][2] = values[17];

        parameters.push_back(param);
        values.clear();
    }
    file.close();
    return parameters;
}

std::vector<CouplingParameter> readCouplingParametersIso(const std::string &filename)
{
    std::vector<CouplingParameter> parameters;
    std::ifstream file(filename);
    std::string line;

    std::cout.precision(20); // does this line only determine the precision of printing or the precision of the string to double conversion?

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return parameters; // Return empty vector if file can't be opened
    }

    // Skip the first line
    std::getline(file, line);

    while (getline(file, line))
    {

        if (line.empty())
        {
            continue;
        }

        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        CouplingParameter param;
        std::istringstream iss(line);
        std::vector<double> values;
        std::string value;

        // retrieve all the values from the line and push them to 'values'
        while (std::getline(iss, value, ','))
        {
            if (!value.empty())
            { // Check if the string is not empty
                try
                {
                    values.push_back(std::stod(value)); // Convert to double and store
                }
                catch (const std::invalid_argument &e)
                {
                    std::cerr << "Invalid argument for std::stod: " << value << std::endl;
                }
                catch (const std::out_of_range &e)
                {
                    std::cerr << "Value out of range for std::stod: " << value << std::endl;
                    continue; // Handle the error
                }
            }
        }

        param.x = values[0];
        param.y = values[1];
        param.z = values[2];

        param.J_iso = values[3];

        parameters.push_back(param);
        values.clear();
    }
    file.close();
    return parameters;
}

std::vector<CouplingParameter> readCouplingParametersPh(const std::string &filename)
{
    std::vector<CouplingParameter> parameters;
    std::ifstream file(filename);
    std::string line;

    std::cout.precision(20); // does this line only determine the precision of printing or the precision of the string to double conversion?

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return parameters; // Return empty vector if file can't be opened
    }

    // Skip the first line
    std::getline(file, line);

    while (getline(file, line))
    {

        if (line.empty())
        {
            continue;
        }

        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        CouplingParameter param;
        std::istringstream iss(line);
        std::vector<double> values;
        std::string value;

        // retrieve all the values from the line and push them to 'values'
        while (std::getline(iss, value, ','))
        {
            if (!value.empty())
            { // Check if the string is not empty
                try
                {
                    values.push_back(std::stod(value)); // Convert to double and store
                }
                catch (const std::invalid_argument &e)
                {
                    std::cerr << "Invalid argument for std::stod: " << value << std::endl;
                }
                catch (const std::out_of_range &e)
                {
                    std::cerr << "Value out of range for std::stod: " << value << std::endl;
                    continue; // Handle the error
                }
            }
        }

        param.x = values[0];
        param.y = values[1];
        param.z = values[2];

        param.Phi[0][0] = values[3];
        param.Phi[0][1] = values[4];
        param.Phi[0][2] = values[5];
        param.Phi[1][0] = values[6];
        param.Phi[1][1] = values[7];
        param.Phi[1][2] = values[8];
        param.Phi[2][0] = values[9];
        param.Phi[2][1] = values[10];
        param.Phi[2][2] = values[11];

        parameters.push_back(param);
        values.clear();
    }
    file.close();
    return parameters;
}

std::string CouplingParameter::printParameter()
{
    std::ostringstream result;

    result << x_ij << ", "
           << y_ij << ", "
           << z_ij << ", "
           << x_jk << ", "
           << y_jk << ", "
           << z_jk << ", "
           << x_ki << ", "
           << y_ki << ", "
           << z_ki << ", "
           << J_xx << ", "
           << J_xy << ", "
           << J_xz << ", "
           << J_yx << ", "
           << J_yy << ", "
           << J_yz << ", "
           << J_zx << ", "
           << J_zy << ", "
           << J_zz;
    return result.str();
}