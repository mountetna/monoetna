import '../lib/client/scss/application.scss'
import {StorybookMockDecorator} from "../lib/client/jsx/storybook-mocks";

export const parameters = {
  actions: { argTypesRegex: "^on[A-Z].*" },
  controls: {
    matchers: {
      color: /(background|color)$/i,
      date: /Date$/,
    },
  },
}

export const decorators = [
  StorybookMockDecorator,
];