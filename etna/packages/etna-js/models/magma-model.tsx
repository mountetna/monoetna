export interface ValidationObject {
  type: string;
  value: string | string[];
}

export interface AttributeObject {
  attribute_name: string;
  attribute_type: string;
  display_name: string;
  description?: string;
  hidden: boolean;
  name: string;
  read_only: boolean;
  restricted: boolean;
  validation: ValidationObject | null;
  link_model_name?: string;
  link_attribute_name?: string;
  link_attribute_type?: string;
  model_name?: string;
}

export interface TemplateObject {
  name: string;
  identifier: string;
  parent: string;
  attributes: {[key: string]: AttributeObject};
}

export interface ModelObject {
  documents: {[key: string]: any};
  revisions: {[key: string]: any};
  template: TemplateObject;
  views: {[key: string]: any};
}

export interface ModelsObject {
  [key: string]: ModelObject;
}

export interface Attribute extends AttributeObject {};

const LINK_TYPES = [ 'collection', 'table', 'child', 'parent', 'link' ];
const LINK_FOREIGN_KEY_TYPES = [ 'parent', 'link' ];
const LINK_COLLECTION_TYPES = [ 'collection', 'table' ];
const LINK_CHILD_TYPES = [ 'collection', 'table', 'child' ];

const FILE_TYPES = [ 'file', 'image', 'file_collection' ];

export class Attribute {
  constructor(attributeDef: AttributeObject) {
    Object.assign(this, attributeDef)
  }

  isType(attributeType: string): boolean;
  isType(attributeTypes: string[]): boolean;
  isType(attributeType: string | string[]): boolean {
    if (typeof attributeType === 'string')
      return this.attribute_type == attributeType;
    else
      return attributeType.includes(this.attribute_type);
  }

  isChildLink(): boolean {
    return LINK_CHILD_TYPES.includes(this.attribute_type);
  }

  isCollectionLink(): boolean {
    return LINK_COLLECTION_TYPES.includes(this.attribute_type);
  }

  isForeignLink(): boolean {
    return LINK_FOREIGN_KEY_TYPES.includes(this.attribute_type);
  }

  isFile(): boolean {
    return this.isType(FILE_TYPES);
  }

}

export interface Model extends TemplateObject {
  attributes: { [key: string]: Attribute }
};

export class Model {
  constructor(modelDef: ModelObject) {
    Object.assign(this, modelDef.template);
    this.attributes = Object.fromEntries(
      Object.entries(modelDef.template.attributes).map(
        ([ attributeName, attributeDef ]) => [
          attributeName, new Attribute(attributeDef)
        ]
      )
    )
  }

  attribute(attributeName: string): Attribute | undefined {
    return this.attributes[ attributeName ];
  }

  selectAttributes(filter: (attribute: Attribute) => boolean): Attribute[] {
    return Object.values(this.attributes).filter(filter);
  }

  hasAttribute(attributeName: string): boolean;
  hasAttribute(filter: (attribute: Attribute) => boolean): boolean;
  hasAttribute(attribute: any): boolean {
    if (typeof attribute === 'string') {
      return attribute in this.attributes;
    }

    return Object.values(this.attributes).some(attribute);
  }

  collects(modelName: string): boolean {
    // ignore self-references
    if (this.name == modelName) return false;

    return this.hasAttribute(
      (a: Attribute) => a.isCollectionLink() && a.link_model_name == modelName
    );
  }

  get isTable(): boolean {
    return this.attributes[ this.parent ]?.link_attribute_type == 'table' 
  }
}

export class Models {
  models: {[key: string]: Model};

  constructor(magmaModels: ModelsObject) {
    this.models = {}
    Object.entries(magmaModels).forEach(
      ([modelName, modelDef]) => {
        this.models[modelName] = new Model(modelDef);
      }
    )
  }

  get size(): number {
    return Object.keys(this.models).length
  }

  get modelNames(): string[] {
    return Object.keys(this.models)
  }

  model(modelName: string): Model | undefined {
    if (modelName in this.models) return this.models[modelName];

    return undefined;
  }

  attribute(modelName: string, attributeName: string): Attribute | undefined {
    if (!(modelName  in this.models)) return undefined;

    return this.models[modelName].attributes[attributeName]
  }

  each( callBack: (modelName:string, model: Model) => void ) {
    Object.entries(this.models).forEach(
      ([modelName, model]: [ modelName: string, model: Model]) => callBack(modelName, model)
    );
  }

  has(modelName: string): boolean {
    return modelName in this.models;
  }
}
