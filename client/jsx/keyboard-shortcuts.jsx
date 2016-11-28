/*
 * Contains the bindings for keyboard shortcuts.
 */
export default class KeyboardShortcuts{

  constructor(){

    document.body.addEventListener('keypress', this.handleKeypress.bind(this));
  }

  handleKeypress(event){

    event = event || window.event;
    //if(event.keyCode == 13 || event.which == 13)
        
    var handled = false;
    var cmd = (event.ctrlKey ? 1 : 0)  |
              (event.altKey ? 2 : 0)   |
              (event.shiftKey ? 4 : 0) |
              (event.metaKey ? 8 : 0);

    /*
    First, handle the key bindings that are independent whether an input
    control is selected or not.
    */
    if(cmd === 1 || cmd === 8 || cmd === 5 || cmd === 12){

      // either CTRL or META key with optional SHIFT.
      switch(event.keyCode){

        case 70: // ctrl + shift + f

          // Find
          break;
        case 71: // ctrl + shift + g

          break;
        case 61:  // FF/Mac '='
        case 107: // FF '+' and '='
        case 187: // Chrome '+'
        case 171: // FF with German keyboard

          // Zoom in or plus
          break;
        case 173: // FF/Mac '-'
        case 109: // FF '-'
        case 189: // Chrome '-'
                
          // Zoom out or subtract
          break;
        case 48: // '0'
        case 96: // '0' on Numpad of Swedish keyboard

          // Reset scale or zoom
          break;
      }
    }
    
    // CTRL or META without shift
    if(cmd === 1 || cmd === 8){
    
      switch(event.keyCode){
    
        case 83: // ctrl + s

          // Save
          break;
        case 14: // ctrl + n

          /*
           * We are using a button to surragate the file input so we may have 
           * a custom browse button.
           */
          document.getElementById('file-selector').click();
          break;
      }
    }
    
    // CTRL+ALT or Option+Command
    if(cmd === 3 || cmd === 10){
    
      switch(event.keyCode){
    
        case 80: //p

          //presentaion mode
          break;
        case 71: //g

          //focus page number dialoge
          break;
      }
    }
  }
}