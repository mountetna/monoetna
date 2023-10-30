import React, { useState, useRef } from "react";



export function inputEventValueToNumber(event: React.ChangeEvent<HTMLInputElement>): number | undefined {
    const eventValue = event.target.value
    return eventValue == "" ? undefined : Number(eventValue)
}