import {listFilesRecursive} from '../file_actions';
import {stubUrl, mockStore} from 'etna-js/spec/helpers';

describe('file_actions', () => {
  describe('listFilesRecursive', () => {
    const bucket_name = 'test-bucket';
    function stubRecursiveListResponse(folderPath, contentsSpec) {
      const joinableFolderPath = (folderPath ? folderPath + '/' : '');
      const subFolders = Object.keys(contentsSpec.folders || {}).map(folder_name => ({
        folder_name,
        folder_path: joinableFolderPath + folder_name,
      }));

      stubUrl({
        verb: 'get',
        path: `/${CONFIG.project_name}/list/${bucket_name}/${folderPath}`,
        response: {
          files: contentsSpec.files.map((file_name) => ({
            file_path: joinableFolderPath + file_name
          })),
          folders: subFolders
        },
        host: 'http://localhost'
      });

      subFolders.forEach(({folder_name, folder_path}) => {
        stubRecursiveListResponse(
          folder_path,
          contentsSpec.folders[folder_name]
        );
      });
    }

    it('flattens out all recursive files belonging to the given root', async () => {
      const folder_name = 'test-folder';
      const store = mockStore({});

      stubRecursiveListResponse(folder_name, {
        files: ['b', 'c', 'a'],
        folders: {
          d: {
            files: ['e', 'g', 'f'],
            folders: {
              h: {
                files: ['i'],
                folders: {}
              }
            }
          },
          j: {
            files: ['k']
          }
        }
      });

      const result = await listFilesRecursive({folder_name, bucket_name})(
        store.dispatch
      );
      expect(result.map(({file_path}) => file_path)).toEqual([
        `${folder_name}/a`,
        `${folder_name}/b`,
        `${folder_name}/c`,
        `${folder_name}/d/e`,
        `${folder_name}/d/f`,
        `${folder_name}/d/g`,
        `${folder_name}/d/h/i`,
        `${folder_name}/j/k`
      ]);
    });
  });
});
