import React, { useState, useRef } from "react";

import { CreateName, Rule } from "../models";



const NameComposeGroup = ({names, rule}: {names: CreateName[], rule: Rule}) => {
    return (
        <div>{`${rule.name} with ${names.length} names`}</div>
    )
}

export default NameComposeGroup;