# AltSplicing Workflow Tutorial

### 1. Run RNA-Seq workflow first

### 2. Login the fat node
```
ssh node2
```

### 3. Clone the repository
```
cd ~/Project/${GROUP}_${DATE}
git clone https://github.com/bioxfu/AltSplicing
cd AltSplicing
```

### 4. Install Spanki (optional)
```
./install_spanki.sh 
```

### 5. Create *workflow.sh* based on the example
```
cp example/workflow.sh workflow.sh

# edit workflow.sh
```

### 6. Initiate the project
```
source init.sh
```

### 7. Run the workflow
```
nohup ./workflow.sh &
```

