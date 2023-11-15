import React from "react";

import NamesToolbar from "../names-toolbar/toolbar";
import AddNamesButton from "../names-toolbar/add-names-button";
import FindAndFilterButton from "../names-toolbar/find-and-filter-button";
import CopyAndReplaceButton from "../names-toolbar/copy-and-replace-button"
import NamesCreateButton from "../names-toolbar/names-create-button";
import DeleteButton from "../names-toolbar/delete-button";



const NamesCreateToolbar = () => {
    return (
        <NamesToolbar
            buttons={[
                <AddNamesButton small={true} />,
                <FindAndFilterButton small={true} />,
                <CopyAndReplaceButton small={true} />,
                <DeleteButton small={true} />,
                <NamesCreateButton small={true} />,
            ]}
        />
    )
};

export default NamesCreateToolbar;