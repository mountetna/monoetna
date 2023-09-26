import pytest

@pytest.fixture
def labor_template():
    return {
        "name": "labor",
        "attributes": {
            "created_at": {
                "name": "created_at",
                "attribute_name": "created_at",
                "display_name": "Created At",
                "hidden": True,
                "validation": None,
                "attribute_type": "date_time"
            },
            "updated_at": {
                "name": "updated_at",
                "attribute_name": "updated_at",
                "display_name": "Updated At",
                "hidden": True,
                "validation": None,
                "attribute_type": "date_time"
            },
            "name": {
                "name": "name",
                "attribute_name": "name",
                "description": "Name for this labor.",
                "display_name": "Name",
                "restricted": False,
                "read_only": False,
                "hidden": False,
                "validation": None,
                "attribute_type": "identifier"
            },
            "country": {
                "name": "country",
                "attribute_name": "country",
                "description": "Place where this labor took place.",
                "restricted": False,
                "read_only": False,
                "hidden": False,
                "validation": None,
                "attribute_type": "string"
            },
            "project": {
                "name": "project",
                "attribute_name": "project",
                "description": "",
                "restricted": False,
                "read_only": False,
                "hidden": False,
                "validation": None,
                "link_model_name": "project",
                "link_attribute_name": "labor",
                "attribute_type": "parent"
            }
        },
        "identifier": "name",
        "parent": "project"
    }

@pytest.fixture
def project_template():
    return {
        "name": "project",
        "attributes": {
            "created_at": {
                "name": "created_at",
                "attribute_name": "created_at",
                "display_name": "Created At",
                "hidden": True,
                "validation": None,
                "attribute_type": "date_time"
            },
            "updated_at": {
                "name": "updated_at",
                "attribute_name": "updated_at",
                "display_name": "Updated At",
                "hidden": True,
                "validation": None,
                "attribute_type": "date_time"
            },
            "name": {
                "name": "name",
                "attribute_name": "name",
                "description": "Name for this project.",
                "display_name": "Name",
                "restricted": False,
                "read_only": False,
                "hidden": False,
                "validation": None,
                "attribute_type": "identifier"
            },
            "labor": {
                "name": "labor",
                "attribute_name": "labor",
                "description": "The Twelve Labors.",
                "restricted": False,
                "read_only": False,
                "hidden": False,
                "validation": None,
                "link_model_name": "labor",
                "link_attribute_name": "project",
                "attribute_type": "collection"
            }
        },
        "identifier": "name",
    }
