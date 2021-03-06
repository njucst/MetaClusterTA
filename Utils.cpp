#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>
#include <getopt.h>

#include "Utils.h"

using namespace std;

struct Option
{
	string long_name;
	string short_name;
	void *pointer;
	OptionType type;
	string description;
	string default_value;
	string value;

	Option(const string &long_name, const string &short_name,
			void *pointer, OptionType type, const string &description, 
			const string &default_value)
	{
		if (short_name.size() > 1 || long_name.size() < 2)
		{
			fprintf(stderr, "invalid parameter: %s, %s\n", 
					long_name.c_str(), short_name.c_str());
			throw exception();
		}

		this->long_name = long_name;
		this->short_name = short_name;
		this->pointer = pointer;
		this->type = type;
		this->description = description;
		this->default_value = default_value;
	}

	void Process()
	{
		if (value == "")
			value = default_value;

		if (value != "")
		{
			switch (type)
			{
				case OptionBool:
					if (value != "" || value == "false")
						*(bool *)pointer = true;
					break;

				case OptionInt:
					*(int *)pointer = atoi(value.c_str());
					break;

				case OptionDouble:
					*(double *)pointer = atof(value.c_str());
					break;

				case OptionString:
					*(string *)pointer = value;
			}
		}
	}

	string ToString()
	{
		string s;

		if (short_name != "")
			s += "  -" + short_name + ", ";
		else
			s += "	  ";

		s += "--" + long_name;
		if (type != OptionBool)
			s += " arg";

		if (default_value != "")
			s += " (=" + default_value + ")";

		while (s.size() < 30)
			s += " ";

		s += " " + description;

		return s;
	}
};

static const int FastaMaxLine = 80;

static map<string, void *> parameters;
static map<string, ParameterType> types;
static vector<Option> options;

void AddOption(const string &long_name, const string &short_name,
		void *pointer, OptionType type, const string &description, const string &default_value)
{
	for (unsigned i = 0; i < options.size(); ++i)
	{
		if (options[i].long_name == long_name
				|| (short_name != "" && options[i].short_name == short_name))
		{
			fprintf(stderr, "parameter already exist: %s, %s\n",
					long_name.c_str(), short_name.c_str());
			throw exception();
		}
	}

	options.push_back(Option(long_name, short_name, pointer, type, description, default_value));
}

void AddOption(const std::string &long_name, const std::string &short_name,
		bool &bool_option, const std::string &description)
{
	AddOption(long_name, short_name, &bool_option, OptionBool, description, "");
}

void AddOption(const std::string &long_name, const std::string &short_name,
		int &int_option, const std::string &description)
{
	stringstream ss;
	ss << int_option;
	AddOption(long_name, short_name, &int_option, OptionInt, description, ss.str());
}

void AddOption(const std::string &long_name, const std::string &short_name,
		double &double_option, const std::string &description)
{
	stringstream ss;
	ss << double_option;
	AddOption(long_name, short_name, &double_option, OptionDouble, description, ss.str());
}

void AddOption(const std::string &long_name, const std::string &short_name,
		std::string &string_option, const std::string &description)
{
	AddOption(long_name, short_name, &string_option, OptionString, description, string_option);
}


void ProcessOptions(int &argc, char *argv[])
{
	struct option long_options[options.size()+1];
	string short_options;
	for (unsigned i = 0; i < options.size(); ++i)
	{
		long_options[i].name = options[i].long_name.c_str();
		long_options[i].flag = 0;
		long_options[i].val = 
			(options[i].short_name != "" ? options[i].short_name[0] : 0);

		if (options[i].type == OptionBool)
		{
			if (options[i].short_name != "")
				short_options += options[i].short_name;
			long_options[i].has_arg = no_argument;
		}
		else
		{
			if (options[i].short_name != "")
				short_options += options[i].short_name + ":";
			long_options[i].has_arg = required_argument;
		}
	}
	long_options[options.size()].name = 0;
	long_options[options.size()].has_arg = 0;
	long_options[options.size()].flag = 0;
	long_options[options.size()].val = 0;

	while (true)
	{
		int index = -1;
		int ch = getopt_long(argc, argv, short_options.c_str(), long_options, &index);

//		if (opterr != 0)
//		{
//			throw exception();
//		}

		if (ch == -1)
			break;

		if (ch == '?')
			throw exception();

		if (ch != 0)
		{
			string s;
			s += ch;

			for (unsigned i = 0; i < options.size(); ++i)
			{
				if (options[i].short_name == s)
				{
					index = i;
					break;
				}
			}
		}

		if (options[index].type == OptionBool)
			options[index].value = "true";
		else
			options[index].value = optarg;
	}

	int index = 1;
	for (int i = optind; i < argc; ++i)
	{
		argv[index++] = argv[i];
	}
	argc = index;

	for (unsigned i = 0; i < options.size(); ++i)
	{
		options[i].Process();
	}
}

string OptionDescriptions()
{
	stringstream ss;
	for (unsigned i = 0; i < options.size(); ++i)
	{
		ss << options[i].ToString() << endl;
	}

	return ss.str();
}


void AddParameter(const char *name, void *pointer, ParameterType type)
{
	parameters[name] = pointer;
	types[name] = type;
}

void ProcessParameters(int &argc, char *argv[])
{
	int k = 0;
	for (int i = 0; i < argc; ++i)
	{
		if (argv[i][0] == '-' && argv[i][1] == '-')
		{
			char *name = argv[i] + 2;
			if (parameters.find(name) != parameters.end())
			{
				switch (types[name])
				{
					case SIMPLE:
						*(bool *)parameters[name] = true;
						break;

					case INTEGER:
						*(int *)parameters[name] = atoi(argv[++i]);
						break;

					case FLOAT:
						*(double *)parameters[name] = atof(argv[++i]);
						break;

					case STRING:
						strcpy((char*)parameters[name], argv[++i]);
				}
			}
		}
		else
		{
			argv[k++] = argv[i];
		}
	}
	argv[k] = NULL;
	argc = k;
}

long long getFileLine(string file)
{
	string command = (string("wc ")+file+" -l > "+file+".tmp");
	system(command.c_str());
	long long ans = 0;
	ifstream ifs((file+".tmp").c_str());
	ifs >> ans;
	ifs.close();
	system((string("rm ")+file+".tmp").c_str());
	return ans;
}
