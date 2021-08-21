import { JSONEditor, defaults } from '@json-editor/json-editor';

JSONEditor.defaults.options.collapsed = true;

export const editConfig = (schema, config) => {
  console.log("Creating editor");
  // cleanup existing editor
  // removeEditor();
 
  // create an editor DOM element
  const editorRoot = document.createElement('div');
  editorRoot.id = 'config-editor';
  document.body.appendChild(editorRoot);
  
  const editor = new JSONEditor( editorRoot, {
    schema, startval: config,
    display_required_only: true,
    theme: 'barebones',
    iconlib: 'fontawesome5',
    remove_button_labels: true,
    compact: true,
    array_controls_top: true,
    disable_array_add: false,
    disable_array_delete: false,
    disable_array_delete_all_rows: true,
    disable_array_delete_last_row: true,
    disable_array_reorder: true,
    disable_edit_json: true,
    disable_properties: false
  });
}
