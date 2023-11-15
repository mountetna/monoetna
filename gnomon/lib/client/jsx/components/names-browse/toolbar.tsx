import React from "react";

import ExportButton from "../names-toolbar/export-button";
import NamesToolbar from "../names-toolbar/toolbar";



const NamesBrowseToolbar = ({ exportData, exportButtonText }: {
    exportData: Array<any>,
    exportButtonText: string,
}) => {
    return (
        <NamesToolbar
            buttons={[
                <ExportButton
                    small={true}
                    data={exportData}
                    buttonText={exportButtonText}
                />
            ]}
        />
    )
};

export default NamesBrowseToolbar;