## Incorporating code from outside repositories

### Examples:

1. Adding an external repository
    ```
    git subtree add --prefix code/Strain-Level_Metagenome_Analysis https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis.git master --squash
    ```

2. Updating an external repository
    ```
    git subtree pull --prefix code/Strain-Level_Metagenome_Analysis https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis.git master --squash
    ```

3. List subtrees
    ```
    git subtree list --resolve
    ```