import * as ReactDOM from "react-dom";
import {Provider} from "react-redux";
import * as React from "react";
import 'regenerator-runtime';
import {VulcanStore} from "./vulcan_store";
import {VulcanProvider} from "./contexts/vulcan_context";

const store = VulcanStore();

function DevApp() {
  return <div>
    <VulcanProvider>
      <div id='ui-container'>
        Helllo
      </div>
    </VulcanProvider>
  </div>
}

window.onload = () => {
  ReactDOM.render(
    <Provider store={store}>
      <DevApp />
    </Provider>,

    document.getElementById("main")
  );
}