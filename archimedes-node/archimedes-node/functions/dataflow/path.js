// This is compiled into JS by the npm script, ``npm run initialize``
//   so it can be imported and used by `index.js`.
import * as _path from 'path';
import * as _fs from 'fs';
import { access as _access, readFile as _readFile, writeFile as _writeFile } from 'fs/promises';
import * as _util from 'util';
import * as _csv from 'csv';
export const inputPath = async (name, inputsEnv = process.env, inputsDir = null) => {
    let path;
    if (null == inputsDir) {
        inputsDir = inputsEnv['INPUTS_DIR'];
    }
    path = _path.normalize(_path.join(inputsDir || '', name));
    return _access(path, _fs.constants.F_OK)
        .then(() => {
        return path;
    })
        .catch((err) => {
        return Promise.reject(`Input key ${name} does not exist for this cell`);
    });
};
export const inputVar = async (name, inputsEnv = process.env, inputsDir = null) => {
    return inputPath(name, inputsEnv, inputsDir).then((path) => {
        return _readFile(path, 'utf8');
    });
};
export const inputTsv = async (name, inputsEnv = process.env, inputsDir = null) => {
    return inputPath(name, inputsEnv, inputsDir)
        .then((path) => {
        return _readFile(path, 'utf8');
    })
        .then((data) => {
        // @ts-ignore: 2554; overloaded function call ... does not match signature
        return _util.promisify(_csv.parse)(data, { delimiter: '\t' });
    });
};
export const inputJson = async (name, inputsEnv = process.env, inputsDir = null) => {
    return inputPath(name, inputsEnv, inputsDir)
        .then((path) => {
        return _readFile(path, 'utf8');
    })
        .then((dataString) => JSON.parse(dataString));
};
export const outputPath = async (name, outputsEnv = process.env, outputsDir = null) => {
    let path;
    if (null == outputsDir) {
        outputsDir = outputsEnv['OUTPUTS_DIR'];
    }
    path = _path.normalize(_path.join(outputsDir || '', name));
    if (outputsEnv['ENFORCE_OUTPUTS_EXIST']) {
        return _access(path, _fs.constants.F_OK)
            .then(() => {
            return path;
        })
            .catch((err) => {
            return Promise.reject(`Output key ${name} does not exist for this cell`);
        });
    }
    else {
        return Promise.resolve(path);
    }
};
export const outputJson = async (data, name, outputsEnv = process.env, outputsDir = null) => {
    return outputPath(name, outputsEnv, outputsDir).then((path) => {
        return _writeFile(path, JSON.stringify(data), { encoding: 'utf8' });
    });
};
// def input_bool(name, inputs_env=_os.environ, inputs_dir=None):
// def str2bool(str):
//     return str.lower() in ("yes", "true", "t", "1")
// return  str2bool(input_var(name, inputs_env, inputs_dir))
// def output_tsv(data, name, outputs_env=_os.environ, outputs_dir=None):
// return data.to_csv( output_path(name, outputs_env, outputs_dir) , sep='\t')
// def output_var(data, name, outputs_env=_os.environ, outputs_dir=None):
// with open(output_path(name, outputs_env, outputs_dir), 'w') as output_file:
//     output_file.write(data)
