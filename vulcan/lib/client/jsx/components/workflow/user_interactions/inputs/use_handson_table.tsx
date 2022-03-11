import React from 'react';

// import editors
import {
  registerEditor,
  BaseEditor,
  HandsontableEditor,
  TextEditor
} from 'handsontable/editors';

// import renderers
import {
  registerRenderer,
  baseRenderer,
  textRenderer
} from 'handsontable/renderers';

// import cell types
import {registerCellType, TextCellType} from 'handsontable/cellTypes';

// import plugins
import {
  ContextMenu,
  CopyPaste,
  DragToScroll,
  Formulas,
  UndoRedo,
  registerPlugin
} from 'handsontable/plugins';

// import translations
import {registerLanguageDictionary, enUS} from 'handsontable/i18n';

export default function useHandsonTable() {
  // register individual translations
  registerLanguageDictionary(enUS);

  // register individual editors
  registerEditor(BaseEditor);
  registerEditor(HandsontableEditor);
  registerEditor(TextEditor);

  // register individual renderers
  registerRenderer(baseRenderer);
  registerRenderer(textRenderer);

  // register individual cell types
  registerCellType(TextCellType);

  // register individual plugins
  registerPlugin(ContextMenu);
  registerPlugin(CopyPaste);
  registerPlugin(DragToScroll);
  registerPlugin(Formulas);
  registerPlugin(UndoRedo);
}
